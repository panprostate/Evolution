# Define continuous, unbounded Cumulative Intensity Functions Lambda(u)
# u = (year + T) / T. Therefore u=0 is onset, u=1 is diagnosis.

Lambda_uniform <- function(u, ...) {
    u
}

# Accelerating: Power-law > 1 (e.g., shape = 2 means linearly increasing rate)
Lambda_accelerating <- function(u, shape = 2, ...) {
    u^shape
}

# Gradual: Decelerating power-law < 1 (e.g., shape = 0.6 means high initial rate that gradually slows)
# (Replacing pbeta, which cannot exceed 1)
Lambda_gradual <- function(u, shape = 0.6, ...) {
    u^shape
}

# Punctuated 
Lambda_punctuated <- function(u, n_bursts = 3, burst_sharpness = 10, ...) {
    # n_bursts: Number of recurrent bursts during the subclonal period [0, T].
    # Using an integer ensures the model perfectly hits Lambda(1) = 1 at diagnosis.
    
    # Calculate which burst cycle we are currently in
    cycle_idx <- floor(u * n_bursts)
    
    # Calculate the fractional progress through the current cycle (0 to 1)
    cycle_progress <- (u * n_bursts) - cycle_idx
    
    # Map the linear progress to a sharp step using a symmetric Beta CDF.
    # High burst_sharpness (e.g., 10) creates long flat periods and sudden vertical jumps.
    step_progress <- pbeta(cycle_progress, shape1 = burst_sharpness, shape2 = burst_sharpness)
    
    # Combine completed cycles with the current step, and scale back down
    res <- (cycle_idx + step_progress) / n_bursts
    
    return(res)
}

# Master function to sample ALL events (subclonal + future) in one continuous process
sample_all_event_times <- function(
    subclonal_events,
    time_subclonal,
    future_horizon,
    model = c("uniform", "gradual", "punctuated", "accelerating"),
    model_params = list(),
    grid_size = 5000
) {
    model <- match.arg(model)
    if (subclonal_events == 0) return(numeric(0))
    
    # Calculate maximum normalized time (u_max)
    u_max <- (time_subclonal + future_horizon) / time_subclonal
    
    # Select the appropriate Lambda function
    lambda_fun <- switch(
        model,
        "uniform"      = Lambda_uniform,
        "gradual"      = Lambda_gradual,
        "punctuated"   = Lambda_punctuated,
        "accelerating" = Lambda_accelerating
    )
    
    # Build grid for inverse transform sampling
    u_grid <- seq(0, u_max * 1.1, length.out = grid_size)
    Lambda_grid <- do.call(lambda_fun, c(list(u = u_grid), model_params))
    
    # Ensure strict monotonicity for approximation
    Lambda_grid <- cummax(Lambda_grid)
    
    # -------------------------------------------------------------
    # Phase 1: Subclonal Events [-T, 0] -> u in [0, 1]
    # We strictly enforce exactly `subclonal_events` in this window
    # -------------------------------------------------------------
    z_sub <- runif(subclonal_events, min = 0, max = 1)
    
    # -------------------------------------------------------------
    # Phase 2: Future Events (0, future_horizon] -> u in (1, u_max]
    # Expected future events based on the continuous rate 
    # -------------------------------------------------------------
    lambda_at_umax <- Lambda_grid[which.min(abs(u_grid - u_max))]
    expected_future <- subclonal_events * (lambda_at_umax - 1)
    
    n_future <- rpois(1, lambda = expected_future)
    if (n_future > 0) {
        z_fut <- runif(n_future, min = 1, max = lambda_at_umax)
    } else {
        z_fut <- numeric(0)
    }
    
    # Combine all events in Lambda space
    z_all <- c(z_sub, z_fut)
    
    # Inverse transform: map z values back to u space
    u_samples <- approx(
        x = Lambda_grid,
        y = u_grid,
        xout = z_all,
        method = "linear",
        ties = "ordered"
    )$y
    
    # Map u back to calendar years: year = u * T - T
    years_samples <- u_samples * time_subclonal - time_subclonal
    
    return(sort(years_samples))
}

simulate_events_monotonic <- function(
    clonal_events,
    subclonal_events,
    time_subclonal,
    years = seq(-5, 15, by = 0.1),
    future_horizon = 15,
    model = c("uniform", "gradual", "punctuated", "accelerating"),
    model_params = list()
) {
    model <- match.arg(model)

    event_times <- sample_all_event_times(
        subclonal_events = subclonal_events,
        time_subclonal = time_subclonal,
        future_horizon = future_horizon,
        model = model,
        model_params = model_params
    )

    total_events <- sapply(years, function(t) {
        clonal_events + sum(event_times <= t)
    })

    data.frame(
        year = years,
        total_events = total_events
    )
}


simulate_pga_uniform <- function(
    pga_diagnosis, 
    yearly_rate, 
    years = seq(-5, 15, by = 0.1), 
    # Scale variance relative to the mean rate (e.g., 10% of the rate)
    # This prevents slow-evolving tumours from having chaotic jumps
    variance_multiplier = 0.10 
) {
    if (pga_diagnosis > 1) { pga_diagnosis <- pga_diagnosis / 100 }
    if (yearly_rate > 1) { yearly_rate <- yearly_rate / 100 }
    
    pga_path <- numeric(length(years))
    idx_t0 <- which.min(abs(years - 0))
    pga_path[idx_t0] <- pga_diagnosis
    
    # 1. FORWARD SIMULATION
    if (idx_t0 < length(years)) {
        for (j in (idx_t0 + 1):length(years)) {
            dt_step <- years[j] - years[j-1]
            mu <- yearly_rate * dt_step
            
            # Dynamic variance: scales naturally with the expected rate
            variance <- (mu * variance_multiplier) * dt_step 
            
            if (mu > 0 && variance > 0) {
                shape <- (mu^2) / variance
                scale <- variance / mu
                increment <- rgamma(1, shape = shape, scale = scale)
            } else {
                increment <- 0
            }
            pga_path[j] <- min(pga_path[j-1] + increment, 1.0)
        }
    }
    
    # 2. BACKWARD SIMULATION
    if (idx_t0 > 1) {
        for (j in (idx_t0 - 1):1) {
            dt_step <- years[j+1] - years[j] 
            mu <- yearly_rate * dt_step
            variance <- (mu * variance_multiplier) * dt_step 
            
            if (mu > 0 && variance > 0) {
                shape <- (mu^2) / variance
                scale <- variance / mu
                increment <- rgamma(1, shape = shape, scale = scale)
            } else {
                increment <- 0
            }
            pga_path[j] <- max(pga_path[j+1] - increment, 0.0)
        }
    }
    
    data.frame(year = years, pga = pga_path)
}
