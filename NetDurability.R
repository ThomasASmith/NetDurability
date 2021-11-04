# Calling script for analyses of PMI net durability data, via Markov model


### SET OPTIONS

# codeMissing indicates how intervals with missing data at the end are treated, Options are:
# 'DXX' code only destroyed or repurposed nets as Attrited
# 'DRX' also code removed nets as Attrited
# 'DRA' also code nets that were absent as Attrited.
# 'DXA' code absent as Attrited, removed as not Attrited .
codeMissing= 'DXA'

# HoleCutoff determines which nets count as damaged
# Nets with pHI <= HoleCutoff are scored as intact
HoleCutoff=20 # 20 for the reference simulation

# To differently code missing data on holes, set value of repairs_recoded to 'intact' or 'holed'
repairs_recoded = 'holed'

# For analysis of a specific net type, change required_net_type to another type in the vector of net types
net_types <- c("Dawaplus 2.0","DuraNet","Interceptor","NetProtect","Olyset","PermaNet 2.0","PermaNet 3.0")
required_net_type='All'

# Number of samples to draw from the posterior densities of MCMC runs
nsamples_posterior =1000 #set to at least 1000 for definitive runs (use smaller values for testing to reduce compute time)

### STANDARD SETTINGS (DO NOT CHANGE THESE)

# Usewts corresponds to reweighting the categories of damage transitions by use categories
reweighting = 'UseWts'

categories <- c("New","Used last night","Not used last night","Not in use","Attrited")
extended_categories <- c("New","Used last night","Not used last night","Not in use","Attrited","Removed","Absent")
hole_categories <- c('New','Undamaged-Not in use','Undamaged-In use','Damaged-In use','Damaged-Not in use','Attrited')
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# define colour palettes for the different stratifications of nets (noting that the 1st and last categories are the same)
modified <- c("#FCFFA4FF","#4DC36BFF","#35608DFF","#CC6666","#000004FF")
modifiedinferno <- c("#FCFFA4FF","#FCAE12FF","#E75E2EFF","#AC2F5EFF","#63156EFF","#000004FF")
# colourBlind safe colors, slightly more blue pink/purple #C879C8 to replace pink#CC79A7
mp_recomCol = c("#D55E00", "#009E73", "#0072A7","#C879C8")

### LOAD ESSENTIAL FUNCTIONS
source("NetDurability_functions.R")
set_graphics_options()

### READ DATA INTO R FROM .csv FILES
survey_data = read.csv('net_status.csv')
all_transitions = read.csv('transitions.csv')
all_interval_estimates = data.frame(Variable=character(0),X2.5.=numeric(0),X50.=numeric(0),X97.5.=numeric(0),
                                    model=character(0),net_type=character(0),repair_recoding=character(0))

transitions <- select_records(all_transitions)
hole_transitions <- get_hole_transitions(transitions=transitions,hc=HoleCutoff)
# Save dataset for fitting of ODE model
write.csv(hole_transitions,file=paste0('hole_transitions',required_options(),'.csv'))

######### RUN THE REQUIRED ANALYSES FROM AMONG THE FOLLOWING:
# Dataset for analysis of whether initial hole index predicts transitions
PHItable <- get_PHItable(transitions=transitions)

## ANALYSIS OF CHANGES IN PHI
# The plot of the results corresponding to Figure 3 is returned as 'PHIScatter.png'
plotPHIscatter()

## FLOW DIAGRAMS
# Figures 4 is returned as 'FlowAgeDependentUseCategoriesAllDXAholedUseWts.png'
# Figures 5 is returned as 'ConstantHoleCategoriesAllDXAholedUseWts.png'
# note that this generates warnings that can be ignored: do not attempt to use geom_flow()
flowdiagrams()

## FITTING OF STAN MODEL
# This can be computationally expensive (depending on the number of iterations, which is specified by parameter iter)
stan_output = ODEFitting(chains = 4,        # number of Markov chains
           warmup = 500,    # number of warmup iterations per chain
           iter = 4000,    # total number of iterations per chain
           thin = 1)        # interval used for thinning outputs (to reduce size of output file))
# The samples from the posterior are stored as a .csv file with default filename 'samples_posteriorAllDXAholedUseWts.csv'
# Convergence diagnostics etc. can be obtained from stan_output$fit_posterior

## COMPUTATION OF INTERVAL ESTIMATES FOR THE ODE PARAMETERS AND DERIVED QUANTITIES
# This requires considerable compute time for nsamples_posterior = 1000
# The plot of the results corresponding to Figure 6 is returned as 'Predictions_AllDXAholedUseWts.png'
interval_estimates = ODEinterval_estimates()

## ANALYSIS OF REPORTED DAMAGE
# The plot of the results corresponding to Figure 7 is returned as 'AnalysisOfRecalls.png'
source("ReportedDamageInAttritedNets.R")
ReportedDamageByAge=getReportedDamageByAge()
ReportedDamageByPreviousDamage = getReportedDamageByPreviousDamage(reported=ReportedDamageByAge$analysis)
ModelForTypeOfAttrition = getModelForTypeOfAttrition(input=ReportedDamageByPreviousDamage)



