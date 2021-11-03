# Function to estimate parameters of ODE model for net durability using rstan
# Written by Adrian Denz Feb 2020.  Converted to a function by Tom Smith to allow external calls
# specifying the data matrix as an argument and returning samples from the posterior
netdurability_run_fixqp <- function(counts=NULL){

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
  
} else {
  counts <- matrix(as.integer(counts),nrow=5,ncol=7)  
}

###################
#models
netdurability_NODATA <- "./netdurability_NODATA.stan"
netdurability_fixqp <- "./netdurability_fixqp.stan"
netdurability_fixq <- "./netdurability_fixq.stan"
netdurability_transformed <- "./netdurability_transformed.stan"
netdurability <- "./netdurability.stan"
###################

#posterior
# data
data_stan <- list(
  # fixed parameters should be assigned here 
  q = 0,
  p0 = 0.602,
  counts = counts,
  t = 0.5,
  priorsigma = 2
)

# fit
fit_posterior<-stan(
  file = netdurability_fixqp,  # Stan program
  data = data_stan,    # named list of data
  chains = 4,          # number of Markov chains
  warmup = 20000,      # number of warmup iterations per chain 
  iter = 1120000,      # total number of iterations per chain
  thin =1000,           # interval used for thinning outputs (to reduce size of output file)
  refresh = 1000,             # progress shown
  control = list(adapt_delta = 0.8)
)
rm(list="data_stan")

pairs(fit_posterior, pars=c('a1', 'a2','a3','a4','h1','h3','u1','u2','v3','v4'))
samples_posterior <- as.data.frame(extract(fit_posterior, 
            pars=c('a1', 'a2','a3','a4','h1','h3','u1','u2','v3','v4')))
#samples_rat <- as.data.frame(extract(fit_posterior, pars=c('rat1', 'rat2')))
setwd('../')
return(samples_posterior)
}

