// Stan model for EACoMoPP semifield experiments, control 
data{
int<lower=0> counts[5,7]; // number of individual nets, rows are S1,...,S4 and columns are S1,...,S4,A,S6,S7
real<lower=0> t; // time between surveys
real<lower=0> priorsigma; // sigma of the lognormal prior on the rates
}
parameters{
  real<lower = 0, upper = 1> p0; // probability that hole status is measured
  real<lower = 0, upper = 1> q; // probability that new net is readily used, i.e. instantly moves to S3, otherwise it moves instantly to S1
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
  matrix<lower = 0, upper=1>[5,7] PP; // probabilities corresponding to data, including transition from new nets and nets with only use measure but not hole measure; BE CAREFUL: rows are for initial states N,S1,...,S4
  // define and compute
  Q = [[-(h1+u1+a1), h1, u1, 0, a1],[0, -(u2+a2), 0, u2, a2],[v3, 0, -(v3+h3+a3), h3, a3],[0, v4, 0, -(v4+a4), a4],[0, 0, 0, 0, 0]];
  P = matrix_exp(t*Q);
  //
  PP[1,1:4] = p0*((1-q)*P[1,1:4] + q*P[3,1:4]); // new nets for which hole status is measured
  PP[1,5] = (1-q)*P[1,5] + q*P[3,5]; // new nets for which hole status is measured
  PP[1,6] = (1-p0)*(((1-q)*P[1,1] + q*P[3,1])+((1-q)*P[1,2] + q*P[3,2])); // new nets for which hole status is measured
  PP[1,7] = (1-p0)*(((1-q)*P[1,3] + q*P[3,3])+((1-q)*P[1,4] + q*P[3,4])); // new nets for which hole status is measured
  //
  PP[2:5,1:4] = p0*P[1:4,1:4]; // non-new nets for which hole status is measured 
  PP[2:5,5] = P[1:4,5]; // non-new nets ending up in A
  PP[2:5,6] = (1-p0)*(P[1:4,1] + P[1:4,2]); // non-new nets for which hole status is measured and that will be used
  PP[2:5,7] = (1-p0)*(P[1:4,3] + P[1:4,4]); // non-new nets for which hole status is measured and that will not be used
}
model{
  // declare
  vector[7] prob; // probabilities
  // priors on rates 
  target += uniform_lpdf(p0 | 0,1);
  target += uniform_lpdf(q | 0,1);
  target += uniform_lpdf(a1 | 0,1);
  target += lognormal_lpdf(a2 | 0, priorsigma);
  target += lognormal_lpdf(a3 | 0, priorsigma);
  target += lognormal_lpdf(a4 | 0, priorsigma);
  target += uniform_lpdf(h1 | 0,1);
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

