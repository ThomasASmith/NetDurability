setwd("~/Documents/R/durabilityanalysis/MarkovmodelStan")
library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# load data
source("./net_load.R")

# Ad-hoc data correction, to be rethought later
# nets that get rid of their holes are excluded
counts[3,1] <- 0;
counts[3,3] <- 0;
counts[5,1] <- 0;
counts[5,3] <- 0;




#posterior
  # data
  data_stan <- list(
    counts = counts,
    t = 0.5,
    priorsigma = 2
  )
  
  # fit
  fit_posterior<-stan(
    file = "./netdurability_uniforma1h1.stan",  # Stan program
    data = data_stan,    # named list of data
    chains = 4,             # number of Markov chains
    warmup = 3000,          # number of warmup iterations per chain
    iter = 6000,            # total number of iterations per chain
    refresh = 0,             # no progress shown
    control = list(adapt_delta = 0.8)
  )



pairs(fit_posterior, pars=c('p0','q','a1','a2','a3','a4','h1','h3', 'u1','u2', 'v3','v4'))

netdurability_uniforma1h1 <- rstan::extract(fit_posterior)
hist(netdurability_uniforma1h1$q, breaks = 100)
hist(netdurability_uniforma1h1$a1, breaks = 100)
hist(netdurability_uniforma1h1$h1, breaks = 100)

write.csv(netdurability_uniforma1h1,file="./netdurability_uniforma1h1_posteriorsample.csv", row.names = FALSE)




#prior
# data i.e. NO-DATA
data_stan <- list(
  t = 0.5,
  priorsigma = 2
)

# fit
fit_prior<-stan(
  file = "./netdurability_NODATA.stan",  # Stan program
  data = data_stan,    # named list of data
  chains = 1,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  refresh = 0,             # no progress shown
  control = list(adapt_delta = 0.8)
)
rm(list="data_stan")

sample_prior <- extract(fit_prior, pars=c('P', 'PP'));

