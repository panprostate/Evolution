list_drivers <- function(trajectories_fps) {
    ms_list = lapply(trajectories_fps, function(x) read.table(x, header = T))
    driver_list = lapply(ms_list, function(ms) {
        unique(paste(ms[,"CNA"], ms[,"ID"], sep='_'))
    })
    return(driver_list)
}

get_trajectory_drivers <- function(trajectories_fps, ppcg_trajectory_fps) {
    ms_list = lapply(trajectories_fps, function(x) read.table(x, header = T))

    ms_list = lapply(ms_list, function(ms) {
        cbind(ms, `cna_id`=apply(ms, 1, function(ms_row) { paste(ms_row["CNA"], ms_row["ID"], sep='_') }))
    })

    num_events = sapply(ms_list, function(ms) {
        length(unique(ms[,"cna_id"]))
    })

    events = sapply(ms_list, function(ms) {
        unique(ms[,"cna_id"])
    }) 

    # Divide by sum of prevalences of events per subset
    ppcg_ms_list = lapply(ppcg_trajectory_fps, function(x) read.table(x, header = T))
    ppcg_ms_list = lapply(ppcg_ms_list, function(ms) {
        cbind(ms, `cna_id`=apply(ms, 1, function(ms_row) { paste(ms_row["CNA"], ms_row["ID"], sep='_') }))
    })
    expected_num_of_events = sapply(ppcg_ms_list, function(ms) {
        # table of cna_ids gets number of samples to have each events_per_sample
        # dividing this by total number of samples gets frequency of each event
        # then we get sum of these frequencies to calculate
        # the expected number of trajectory events per sample
        sum(table(ms$cna_id) / length(unique(ms$Tumour_Name)))
    })
    
    # Get all unique samples across all trajectories to ensure consistent output
    all_samples = unique(unlist(lapply(ms_list, function(ms) unique(ms[,"Tumour_Name"]))))

    total_events_per_sample = lapply(ms_list, function(ms) {
        sample_counts = table(ms[,"Tumour_Name"])
        # Add 0 counts for samples not present in this trajectory
        missing_samples = setdiff(unique(ms[,"Tumour_Name"]), names(sample_counts))
        if (length(missing_samples) > 0) {
            missing_counts = rep(0, length(missing_samples))
            names(missing_counts) = missing_samples
            sample_counts = c(sample_counts, missing_counts)
        }
        sample_counts[unique(ms[,"Tumour_Name"])]  # Return in consistent order
    })

    clonal_events_per_sample = lapply(ms_list, function(ms) {
        clonal_events = ms[ms[,"clonality"] == "clonal",]
        sample_counts = table(clonal_events[,"Tumour_Name"])
        # Add 0 counts for samples not present in clonal events
        missing_samples = setdiff(unique(ms[,"Tumour_Name"]), names(sample_counts))
        if (length(missing_samples) > 0) {
            missing_counts = rep(0, length(missing_samples))
            names(missing_counts) = missing_samples
            sample_counts = c(sample_counts, missing_counts)
        }
        sample_counts[unique(ms[,"Tumour_Name"])]  # Return in consistent order
    })

    subclonal_events_per_sample = lapply(ms_list, function(ms) {
        subclonal_events = ms[ms[,"clonality"] == "subclonal",]
        sample_counts = table(subclonal_events[,"Tumour_Name"])
        # Add 0 counts for samples not present in subclonal events
        missing_samples = setdiff(unique(ms[,"Tumour_Name"]), names(sample_counts))
        if (length(missing_samples) > 0) {
            missing_counts = rep(0, length(missing_samples))
            names(missing_counts) = missing_samples
            sample_counts = c(sample_counts, missing_counts)
        }
        sample_counts[unique(ms[,"Tumour_Name"])]  # Return in consistent order
    })

    poteo_per_sample = lapply(seq_along(total_events_per_sample), function(subset_i) {
        total_events_per_sample[[subset_i]] / num_events[[subset_i]]
    })
    
    poteo_clonal_per_sample = lapply(seq_along(clonal_events_per_sample), function(subset_i) {
        clonal_events_per_sample[[subset_i]] / num_events[[subset_i]]
    })

    poteo_subclonal_per_sample = lapply(seq_along(subclonal_events_per_sample), function(subset_i) {
        subclonal_events_per_sample[[subset_i]] / num_events[[subset_i]]
    })

    poete_per_sample = lapply(seq_along(total_events_per_sample), function(subset_i) {
        total_events_per_sample[[subset_i]] / expected_num_of_events[[subset_i]]
    })

    poete_clonal_per_sample = lapply(seq_along(clonal_events_per_sample), function(subset_i) {
        clonal_events_per_sample[[subset_i]] / expected_num_of_events[[subset_i]]
    })

    poete_subclonal_per_sample = lapply(seq_along(subclonal_events_per_sample), function(subset_i) {
        subclonal_events_per_sample[[subset_i]] / expected_num_of_events[[subset_i]]
    })

    return(
        data.frame(
            Tumour_Name = names(unlist(total_events_per_sample)), 
            total_events = unlist(total_events_per_sample),
            clonal_events = unlist(clonal_events_per_sample),
            subclonal_events = unlist(subclonal_events_per_sample),
            total_poteo = unlist(poteo_per_sample),
            clonal_poteo = unlist(poteo_clonal_per_sample),
            subclonal_poteo = unlist(poteo_subclonal_per_sample), 
            total_poete = unlist(poete_per_sample),
            clonal_poete = unlist(poete_clonal_per_sample),
            subclonal_poete = unlist(poete_subclonal_per_sample)
        )
    )
}
