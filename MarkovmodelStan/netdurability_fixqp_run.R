# Function to estimate parameters of ODE model for net durability using rstan
# Written by Adrian Denz Feb 2020.  Converted to a function by Tom Smith to allow external calls
# specifying the data matrix as an argument and returning samples from the posterior
netdurability_run_fixq <- function(counts=NULL){

setwd('../MarkovmodelStan')
library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = FALSE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')

if(is.null(counts)){
  # load data
  # produce matrix of counts
  counts <- matrix(data=0,nrow=5,ncol=7)
  source("./net_load.R")
  
  # Ad-hoc data correction, to be rethought later
  # nets that get rid of their holes are excluded
  counts[3,1] <- 0;
  counts[3,3] <- 0;
  counts[5,1] <- 0;
  counts[5,3] <- 0;
} else {
  counts <- matrix(as.integer(counts),nrow=5,ncol=7)  
}

###################
#models
netdurability_NODATA <- "./netdurability_NODATA.stan"
netdurability_fixq <- "./netdurability_fixq.stan"
netdurability <- "./netdurability.stan"
###################

#posterior
# data
data_stan <- list(
  # q should be assigned in the case of the fixed q model 
  q = 0,
  counts = counts,
  t = 0.5,
  priorsigma = 2
)

# fit
fit_posterior<-stan(
  file = netdurability_fixq,  # Stan program
  data = data_stan,    # named list of data
  chains = 8,          # number of Markov chains
  warmup = 2000,      # number of warmup iterations per chain 
  iter = 22000,      # total number of iterations per chain
  thin =1,           # interval used for thinning outputs (to reduce size of output file)
  refresh = 1000,             # progress shown
  control = list(adapt_delta = 0.8)
)
rm(list="data_stan")

pairs(fit_posterior, pars=c('p0','a1', 'a2','a3','a4','h1','h3','u1','u2','v3','v4'))
samples_posterior <- as.data.frame(extract(fit_posterior, 
            pars=c('p0','a1', 'a2','a3','a4','h1','h3','u1','u2','v3','v4')))
setwd('../')
return(samples_posterior)
}

