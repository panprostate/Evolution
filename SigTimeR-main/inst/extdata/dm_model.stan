data {
  int<lower = 0> S;
  int<lower = 0> N;
  int<lower = 0> Y[N];
  vector<lower=0>[S] alpha;
  matrix[N, S] signatures;
}

parameters {
  simplex[S] weights;
}

model {
  weights ~ dirichlet(alpha);
  Y ~ multinomial(signatures * weights);
}
