# Calling script for analyses of PMI net durability data, via Markov model
# The options to select specific analyses are by default set to FALSE 
#
# select which analyses are required
requireDataDescription = FALSE
requireODEFitting = FALSE
requireODEpostprocessing =FALSE
requireCountryComparisons = FALSE
requireNet_typeComparisons = FALSE
requireCalculationIntervalEstimates = FALSE
requireRecyclingAnalysis = FALSE
# codeMissing indicates how intervals with missing data at the end are treated, Options are:
# 'DXX' code only destroyed or repurposed nets as Attrited
# 'DRX' also code removed nets as Attrited
# 'DRA' also code nets that were absent as Attrited.
# 'DXA' code absent as Attrited, removed as not Attrited .
codeMissing= 'DXA'
# Nets with pHI <= HoleCutoff are scored as intact
HoleCutoff=20 # 20 for the reference simulation
# to generate a dataset for fitting (note there is a different function for each choice about how
# to deal with different missing data on holes, set value of repairs_recoded to 'intact' or 'holed')
repairs_recoded = 'holed'
# Usewts corresponds to reweighting the categories of damage transitions by use categories
categories <- c("New","Used last night","Not used last night","Not in use","Attrited")
extended_categories <- c("New","Used last night","Not used last night","Not in use","Attrited","Removed","Absent")
hole_categories <- c('New','Undamaged-Not in use','Undamaged-In use','Damaged-In use','Damaged-Not in use','Attrited')
reweighting = 'UseWts'
required_net_type='All'
required_country='Kenya'
# All possible net types accross all countries are:
# net_types <- c("Dawaplus 2.0","DuraNet","Interceptor","LifeNet","NetProtect","Olyset","PermaNet 2.0","PermaNet 3.0")
net_types <- c("Dawaplus 2.0","DuraNet","Interceptor","NetProtect","Olyset","PermaNet 2.0","PermaNet 3.0")
#select required countries here
countries <- c('Kenya')
nsamples_posterior =1000 #set to 1000 for definitive runs (use smaller values for testing)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
folderForDurationFiles = ''
pathToSavePlots = paste0(folderForDurationFiles,'temp/')

# define separate palettes for the different stratifications of nets (noting that the 1st and last categories are the same)
modified <- c("#FCFFA4FF","#4DC36BFF","#35608DFF","#CC6666","#000004FF")
modifiedinferno <- c("#FCFFA4FF","#FCAE12FF","#E75E2EFF","#AC2F5EFF","#63156EFF","#000004FF")
# colourBlind safe colors, slightly more blue pink/purple #C879C8 to replace pink#CC79A7
mp_recomCol = c("#D55E00", "#009E73", "#0072A7","#C879C8")
source("PMI_Markov_functions.R")
set_graphics_options()
DataWide <- readPMIdataForLongitudinalAnalysis()
all_interval_estimates = data.frame(Variable=character(0),X2.5.=numeric(0),X50.=numeric(0),X97.5.=numeric(0),
                                    model=character(0),net_type=character(0),country=character(0),repair_recoding=character(0))
for(required_country in countries){
  print(required_country)
  if(requireDataDescription | requireODEFitting) {
    dfw <- select_records(DataWide)
    transitions <- get_transitions(dfw=dfw)
    hole_transitions <- get_hole_transitions(transitions=transitions,hc=HoleCutoff)
    plotPHIscatter()
    plotReportedDamage()
  }
  interval_estimates <- processODE(required_net_type=required_net_type,
                                   required_country=required_country)
  all_interval_estimates <- rbind(all_interval_estimates,interval_estimates)
}

if (requireRecyclingAnalysis) {
  model4 <- RecyclingAnalysis()
  model3 <- SensitivityAnalysis()
}

if (requireNet_typeComparisons){
  cumulated_interval_estimates <-data.frame(matrix(ncol=8,nrow=0, dimnames=list(NULL, c("Variable","X2.5.","X50.","X97.5.","model","net_type","country","repair_recoding"))))
  #  for (required_country in c('All',countries)){
  for (required_country in countries){
    print(required_country)
    for(required_net_type in c('All',net_types)){
      print(required_net_type)
      dfw <- select_records(DataWide)
      transitions <- get_transitions(dfw=dfw)
      hole_transitions <- get_hole_transitions(transitions=transitions,hc=HoleCutoff)
      if(nrow(hole_transitions) > 0){
        if(!is.null(interval_estimates)){
          interval_estimates <- processODE(required_net_type=required_net_type,required_country=required_country)
          all_interval_estimates <- read.csv(file=paste0('interval_estimates',required_options(),'.csv'))
          vnames <- names(interval_estimates)
          all_interval_estimates <- all_interval_estimates[,vnames]
          all_interval_estimates <- rbind(all_interval_estimates,interval_estimates)
          write.csv(all_interval_estimates,file=paste0('interval_estimates',required_options(),'.csv'))
          cumulated_interval_estimates <- rbind(cumulated_interval_estimates,interval_estimates)
        }
      }
    }
  }
  plotComparisons(all_interval_estimates=cumulated_interval_estimates,classifier='net_type',filterval='Kenya')
  cumulated_interval_estimates$grouping=with(cumulated_interval_estimates,paste(country,net_type,repair_recoding))
  tabulate_interval_estimates(cumulated_interval_estimates,outfile='EstimatesByNetType.csv')
}

