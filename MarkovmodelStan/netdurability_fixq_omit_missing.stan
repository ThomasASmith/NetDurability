// Stan model for net durability surveys with q fixed 
data{
int<lower=0> counts[5,5]; // number of individual nets, rows are S1,...,S4 and columns are S1,...,S4,A
real<lower=0> t; // time between surveys
real<lower=0> priorsigma; // sigma of the lognormal prior on the rates
}
parameters{
  real<lower=0> a1; // see diagram
  real<lower=0> a2; // see diagram
  real<lower=0> a3; // see diagram
  real<lower=0> a4; // see diagram
  real<lower=0> h1; // see diagram
  real<lower=0> h3; // see diagram
  real<lower=0> u1; // see diagram
  real<lower=0> u2; // see diagram
  real<lower=0> v3; // see diagram
  real<lower=0> v4; // see diagram
}
transformed parameters{
  // declare
  matrix[5,5] Q; // infinitesimal generator matrix for transitions in markov chain on state space S1,...,S4,A
  matrix<lower = 0, upper=1>[5,5] P; // transition probability function; BE CAREFUL: rows are for states S1,...,S4,A
  matrix<lower = 0, upper=1>[5,5] PP; // probabilities corresponding to data; BE CAREFUL: rows are for initial states N,S1,...,S4
  // define and compute
  Q = [[-(h1+u1+a1), h1, u1, 0, a1],[0, -(u2+a2), 0, u2, a2],[v3, 0, -(v3+h3+a3), h3, a3],[0, v4, 0, -(v4+a4), a4],[0, 0, 0, 0, 0]];
  P = matrix_exp(t*Q);
  //
  PP[1,1:4] = P[1,1:4] ; // new nets for which hole status is measured
  PP[1,5] = P[1,5] ; // new nets for which hole status is measured
  //
  PP[2:5,1:4] = P[1:4,1:4]; // non-new nets for which hole status is measured 
  PP[2:5,5] = P[1:4,5]; // non-new nets ending up in A
}
model{
  // declare
  vector[5] prob; // probabilities
  // priors on rates 
  target += lognormal_lpdf(a1 | 0, priorsigma);
  target += lognormal_lpdf(a2 | 0, priorsigma);
  target += lognormal_lpdf(a3 | 0, priorsigma);
  target += lognormal_lpdf(a4 | 0, priorsigma);
  target += lognormal_lpdf(h1 | 0, priorsigma);
  target += lognormal_lpdf(h3 | 0, priorsigma);
  target += lognormal_lpdf(u1 | 0, priorsigma);
  target += lognormal_lpdf(u2 | 0, priorsigma);
  target += lognormal_lpdf(v3 | 0, priorsigma);
  target += lognormal_lpdf(v4 | 0, priorsigma);
 // data model 
  for (i in 1:5){
    prob = (PP[i,])';
    target += multinomial_lpmf(counts[i,]| prob);
  }
}

