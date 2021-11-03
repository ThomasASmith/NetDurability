# Functions called from by PMI_Markov.R
# Code assembled by Tom Smith for cohort analysis of LLIN durability
# February 2020 - November 2020

######################### READING DATA ##############################
# Read raw data and code the required status variables
read_pooled_llin_data <- function(){
  #suppress warning about missing data
  suppressWarnings(df0 <- read.csv("../pooled_llin.sas7bdat.csv"))
  df0=df0[!is.na(df0$PBrand),] 
  df0=df0[df0$PCountry !='Rwanda' & df0$PCountry !='' & df0$PCountry !='Mozambique' 
          & df0$PCountry !='Angola' & df0$PCountry !='Malawi' & df0$PCountry !='Benin' & df0$PCountry !='Zambia',]
  if(required_country !='All'){
    df0=df0[df0$PCountry == required_country,]
  }
  Country <- droplevels(df0$PCountry)
  # CLASSIFICATION OF STATUS AT TIME OF SURVEY
  Status = as.character(with(df0,ifelse(PLLINUse=='None','Not in use','Not used last night')))
  
  # Values for PReasonAbsent are:
  #[1] ""                          "1"                         "Being Washed"             
  #[4] "Damaged"                   "Destroyed"                 "Destroyed--Burned by fire"
  #[7] "Discarded--Too torn up"    "Discarded-Not killing mos" "Do Not Know"              
  #[10] "Eaten by Rats"             "Lost/Stolen"               "Other"                    
  #[13] "Other Use"                 "Replaced (New/Separate IT" "Sold/Given away"          
  #[16] "Stored"                    "Taken from house/Moved" 
  

  Status[as.character(df0$PReasonAbsent) == 'Destroyed--Burned by fire' | 
           as.character(df0$PReasonAbsent) == 'Discarded--Too torn up' | 
           as.character(df0$PReasonAbsent) == 'Discarded-Not killing mos'] <- 'Attrited'
  Status[as.character(df0$PReasonAbsent) == 'Other' | 
           as.character(df0$PReasonAbsent) == 'Other Use'] <- 'Attrited'
  Status[as.character(df0$PReasonAbsent) == 'Taken from house/Moved' | 
         as.character(df0$PReasonAbsent) == 'Lost/Stolen' | 
         as.character(df0$PReasonAbsent) == 'Sold/Given away'] <- 'Removed'
  # levels(df0$PLocationLLIN)
  # [1] ""                     "189"                  "Do Not Know"          "Folded Up"           
  # [5] "Hanging"              "Hanging - Loose"      "Hanging - Tied"       "Inside - Not Hanging"
  # [9] "Not Folded Up"        "Not Hanging"          "Stored Away" 
  
  # Non missing values of PLocationLLIN also indicate that the net is present
  Status[as.numeric(df0$PLocationLLIN) > 3 & Status=='Attrited'] <- 'Not in use'
  # If used last night, then Status must be in use
  Status[df0$PLastNight==1] <- 'Used last night'
  
  # Assign PHI3_Field as default measure of holed area (Kilian's Proportional Hole Index).  
  # Use lab data if these are all that is available (note that scatter indicates acceptable agreement)
  
  PHI <- with(df0,ifelse(is.na(PHI3_Field),PHI3_Lab,PHI3_Field))
  PHI[Status=='Attrited'] <- NA
  
  # compute indicator of whether non-use or attrition is reported as due to damage
  Damage_related <- as.numeric(df0$PReasonAbsent)
  Damage_related <- with(df0, ifelse(Damage_related > 3 & Damage_related < 7 & Damage_related !=5,'Yes','No'))
  
  library(car)
  # Recode wash frequency as a number
  WashFreq <- as.numeric(as.character((car::recode(df0$PWashFreq,"''=NA;'1 Per Year'=1;
'10 Per Year'=10;'12 Per Year'=12;'13-20 Per Year'=17;
'2 Per Year'=2;'21+ Per Year'=22;'3 Per Year'=3;
'4 Per Year'=4;'5 Per Year'=5;'6 Per Year'=6;
'7 Per Year'=7;'8 Per Year'=8;'Never'=0;'Unknown'=NA"))))
  
  # Reduce categories for bedding
  Bedding <- car::recode(df0$PBedding,"''=NA; 'Bamboo'='Other';'Bed'='Mattress/Bed';'Blankets/Sacks/Tarp/Plast'='Mat etc';
'Brick'='Other';'Frame - Metal'='Frame';'Frame - Wood'='Frame';'Ground'='Other';'Mat'='Mat etc';
'Mattress'='Mattress/Bed';'Other'='Other';'Sofa'='Other';'Unknown'=NA")
  
  #
  suppressWarnings(Chem<-as.numeric(as.character(df0$mg_per_m2)))
  df<- with(df0,data.frame(Country=Country,NetCode=PNetCode,Brand=PBrand,Round=PRound,Status=Status,
                           PHI=PHI,Damage_related=Damage_related,Wash=PWash,WashFreq=WashFreq,
                           Bedding=Bedding,Reason=PReasonAbsent,Location=PLocationLLIN))
  return(df)
}


# Add hole status variable to cross sectional data frame
get_holeStatus = function(df,hc=HoleCutoff) {
  df$HoleStatus <-  with(df,ifelse(PHI<=hc & Status !='Not in use','Undamaged-In use',NA))
  df$HoleStatus <-  with(df,ifelse(PHI > hc & Status !='Not in use','Damaged-In use',HoleStatus.i))
  df$HoleStatus <-  with(df,ifelse(PHI<=hc & Status =='Not in use','Undamaged-Not in use',HoleStatus.i))
  df$HoleStatus <-  with(df,ifelse(PHI > hc & Status =='Not in use','Damaged-In use',HoleStatus.i))
  hole_transitions$longData=ifelse(as.numeric(as.factor(hole_transitions$HoleStatus.f)) %in% c(2,3,4,5),0,1)
return(df)}


# Read and merge-in auxiliary data required for linking records longitudinally
readPMIdataForLongitudinalAnalysis <- function(){
  df <- read_pooled_llin_data()
  df2 <- read.csv("netIDresolution.csv")
  df2 <- df2[,c("Country","Brand","NetCode","NetID")]
  df3 <- merge(df,df2)
#  df3$Status[df3$Status=='Removed']<-'Attrited'
  df3$id <- with(df3,paste0(Country,Brand,NetID))
  dfw <- reshape(data=df3,idvar="id",
                 timevar = "Round",
                 v.names = c("Status","PHI","Damage_related","Wash","WashFreq","Bedding"),
                 direction="wide")
  dfw <- dfw[!is.na(dfw$NetID),]
  # dfw is a data frame with one record per net included in longitudinal studies
return(dfw)}

# Dataset for analysis of transitions, not considering whether holes are present
get_transitions <- function(dfw=dfw){
  # Create long format with records for transitions
  dummy <- factor(0)
  length(dummy) <- 0
  # To avoid coding the follow-up from the last round as attrition, check if Age.f is greater than the maximum
  maxAges <- function(Country=Country){
    # make lookup vector with names
    countries = c("Kenya","Senegal")
    maxAges = c(48,36)
    lookup = match(Country, countries)
  return(maxAges[lookup])} 
  
    last <- maxAges(dfw$Country)
    
  transitions0 = data.frame(Country=character(0),NetID=character(0),Brand=character(0),Age.f=numeric(0),
               Status.i=dummy,PHI.i=numeric(0),Status.f=dummy,PHI.f=numeric(0),Damage_related=dummy,Wash=integer(0),WashFreq=numeric(0),Bedding=dummy)
  if(with(dfw,!exists('Status.48'))) { transitions8 <- transitions0 } else { 
    last <- ifelse(is.na(dfw$Status.48) & last==48, 42, last)  
    transitions8 <- with(dfw,data.frame(Country=dfw$Country,NetID=dfw$NetID,Brand=dfw$Brand,Age.f=48,Status.i=Status.42,PHI.i=PHI.42,Status.f=Status.48,PHI.f=PHI.48,
                                      Damage_related=Damage_related.48,Wash=Wash.48,WashFreq=WashFreq.48,Bedding=Bedding.48))}
  if(with(dfw,!exists('Status.42'))) { transitions7 <- transitions0 } else { 
    last <- ifelse(is.na(dfw$Status.42) & last==42, 36, last) 
    transitions7 <- with(dfw,data.frame(Country=dfw$Country,NetID=dfw$NetID,Brand=dfw$Brand,Age.f=42,Status.i=Status.36,PHI.i=PHI.36,Status.f=Status.42,PHI.f=PHI.42,
                                      Damage_related=Damage_related.42,Wash=Wash.42,WashFreq=WashFreq.42,Bedding=Bedding.42))}
  if(with(dfw,!exists('Status.36'))) { transitions6 <- transitions0 } else { 
    last <- ifelse(is.na(dfw$Status.36) & last==36, 30, last) 
    transitions6 <- with(dfw,data.frame(Country=dfw$Country,NetID=dfw$NetID,Brand=dfw$Brand,Age.f=36,Status.i=Status.30,PHI.i=PHI.30,Status.f=Status.36,PHI.f=PHI.36,
                                      Damage_related=Damage_related.36,Wash=Wash.36,WashFreq=WashFreq.36,Bedding=Bedding.36))}
  if(with(dfw,!exists('Status.30'))) { transitions5 <- transitions0 } else { 
    last <- ifelse(is.na(dfw$Status.30) & last==30, 24, last) 
    transitions5 <- with(dfw,data.frame(Country=dfw$Country,NetID=dfw$NetID,Brand=dfw$Brand,Age.f=30,Status.i=Status.24,PHI.i=PHI.24,Status.f=Status.30,PHI.f=PHI.30,
                                        Damage_related=Damage_related.30,Wash=Wash.30,WashFreq=WashFreq.30,Bedding=Bedding.30))}
  if(with(dfw,!exists('Status.24'))) { transitions4 <- transitions0 } else { 
    last <- ifelse(is.na(dfw$Status.24) & last==24, 18, last) 
    transitions4 <- with(dfw,data.frame(Country=dfw$Country,NetID=dfw$NetID,Brand=dfw$Brand,Age.f=24,Status.i=Status.18,PHI.i=PHI.18,Status.f=Status.24,PHI.f=PHI.24,
                                        Damage_related=Damage_related.24,Wash=Wash.24,WashFreq=WashFreq.24,Bedding=Bedding.24))}
  if(with(dfw,!exists('Status.18'))) { transitions3 <- transitions0 } else { 
    last <- ifelse(is.na(dfw$Status.18) & last==18, 12, last) 
    transitions3 <- with(dfw,data.frame(Country=dfw$Country,NetID=dfw$NetID,Brand=dfw$Brand,Age.f=18,Status.i=Status.12,PHI.i=PHI.12,Status.f=Status.18,PHI.f=PHI.18,
                                        Damage_related=Damage_related.18,Wash=Wash.18,WashFreq=WashFreq.18,Bedding=Bedding.18))}
  if(with(dfw,!exists('Status.12'))) { transitions2 <- transitions0 } else { 
    last <- ifelse(is.na(dfw$Status.12) & last==12, 6, last) 
    transitions2 <- with(dfw, data.frame(Country=dfw$Country,NetID=dfw$NetID,Brand=dfw$Brand,Age.f=12,Status.i=Status.6,PHI.i=PHI.6,Status.f=Status.12,PHI.f=PHI.12,
                                      Damage_related=Damage_related.12,Wash=Wash.12,WashFreq=WashFreq.12,Bedding=Bedding.12))}
  if(with(dfw,!exists('Status.6'))) { transitions1 <- transitions0 } else { 
    last <- ifelse(is.na(dfw$Status.6) & last==6, 0, last) 
    transitions1 <- with(dfw, data.frame(Country=dfw$Country,NetID=dfw$NetID,Brand=dfw$Brand,Age.f=6,Status.i='New',PHI.i=0,Status.f=Status.6,PHI.f=PHI.6,
                                         Damage_related=Damage_related.6,Wash=Wash.6,WashFreq=WashFreq.6,Bedding=Bedding.6))}
  
  transitions <- rbind(transitions1,transitions2,transitions3,transitions4,transitions5,transitions6,transitions7,transitions8)
  transitions$Status.f <- as.character(transitions$Status.f)
  transitions$Status.f[is.na(transitions$Status.f)] <- 'Absent'
  # remove uninformative records obtain a dataset for analysis of effect of initial use on final status
  if(codeMissing=='DRA' | codeMissing=='DXA') {
    # Assign as attrited all nets that disappear before the last survey for the country
    transitions$Status.f <- with(transitions,ifelse(is.na(Status.f) & (Age.f==last+6) & (Age.f < maxAges(Country)), 'Absent', Status.f))    
  }
  transitions$ExtendedStatus.f <- factor(transitions$Status.f,levels=extended_categories)
  transitions$Status.f[transitions$Status.f == 'Removed'] <- ifelse(codeMissing=='DXX' | codeMissing=='DXA', NA, 'Attrited')
  transitions$Status.f[transitions$Status.f == 'Absent'] <- ifelse(codeMissing=='DRX' | codeMissing=='DXX', NA, 'Attrited')
  
  transitions <- transitions[!is.na(transitions$Status.i) & !is.na(transitions$Status.f),]
  transitions <- transitions[transitions$Status.i != 'Attrited',]
  transitions <- transitions[transitions$Status.i != 'Removed',]
  transitions <- transitions[transitions$Status.f != 'Unclear',]
  transitions$Status.i <- droplevels(transitions$Status.i)

  transitions$Status.i <- factor(transitions$Status.i,levels=categories)
  transitions$Status.f <- droplevels(as.factor(transitions$Status.f))
  transitions$Status.f <- factor(transitions$Status.f,levels=categories)
  return(transitions)}

# Dataset for analysis of transitions, including whether holes are present
get_hole_transitions <- function(transitions,hc=0){
  hole_transitions <- with(transitions, transitions[!is.na(PHI.i) | Status.i=='New',])
  hole_transitions$HoleStatus.i <-  with(hole_transitions,ifelse(PHI.i<=hc & Status.i !='Not in use','Undamaged-In use',NA))
  hole_transitions$HoleStatus.i <-  with(hole_transitions,ifelse(PHI.i > hc & Status.i !='Not in use','Damaged-In use',HoleStatus.i))
  hole_transitions$HoleStatus.i <-  with(hole_transitions,ifelse(PHI.i<=hc & Status.i =='Not in use','Undamaged-Not in use',HoleStatus.i))
  hole_transitions$HoleStatus.i <-  with(hole_transitions,ifelse(PHI.i > hc & Status.i =='Not in use','Damaged-Not in use',HoleStatus.i))
  hole_transitions$HoleStatus.i <-  with(hole_transitions,ifelse(Status.i =='New','New',HoleStatus.i))
  hole_transitions$HoleStatus.f <-  with(hole_transitions,ifelse(PHI.f<=hc & Status.f !='Not in use','Undamaged-In use',NA))
  hole_transitions$HoleStatus.f <-  with(hole_transitions,ifelse(PHI.f > hc & Status.f !='Not in use','Damaged-In use',HoleStatus.f))
  hole_transitions$HoleStatus.f <-  with(hole_transitions,ifelse(PHI.f<=hc & Status.f =='Not in use','Undamaged-Not in use',HoleStatus.f))
  hole_transitions$HoleStatus.f <-  with(hole_transitions,ifelse(PHI.f > hc & Status.f =='Not in use','Damaged-Not in use',HoleStatus.f))
  hole_transitions$HoleStatus.f <-  with(hole_transitions,ifelse(is.na(PHI.f) & Status.f !='Not in use','NA-In use',HoleStatus.f))
  hole_transitions$HoleStatus.f <-  with(hole_transitions,ifelse(is.na(PHI.f) & Status.f =='Not in use','NA-Not in use',HoleStatus.f))
  hole_transitions$HoleStatus.f <-  with(hole_transitions,ifelse(Status.f =='Attrited','Attrited',HoleStatus.f))
  hole_transitions$longData=ifelse(as.numeric(as.factor(hole_transitions$HoleStatus.f)) %in% c(2,3,4,5),1,0)
  return(hole_transitions)}

# Option to fill missing values for bedding assuming these are invariant for any net
# (this is not implemented as there appears to be information in the changes in washing with age)
fillMissing <- function(stem,dfw=dfw){
  pval=rep(NA,nrow(dfw))
  for (i in 1:8){
    vname <- paste0('dfw$',stem,'.',6*i)
    if(is.na(eval(parse(text = vname)))) eval(parse(text = vname)) <- pval
    pval <- eval(parse(text = vname))
  }
}

# Dataset for analysis of whether initial hole index predicts transitions 
get_PHItable <- function(transitions=transitions){
  analysis1 <- transitions[!is.na(transitions$PHI.i) 
                           & transitions$Status.i !='Removed'
                           & transitions$Status.i !='New',]
  analysis1$Status.i <- droplevels(analysis1$Status.i)
  PHItable <- with(analysis1,table(Status.i,Status.f))
  return(PHItable)
}

set_graphics_options <- function() {
  # Graphics options for plots of vector control parameterisations
  # Use for plotting of analyses of semi-field experiments and of field studies of durability
  # written by Tom Smith, 2019-2020
  library(xlsx)
  library(gdata)
  library(RCurl)
  library(reshape2)
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(scales)
  
}


# Function for saving plots
savePlot <- function(plot,Plotname,vertical_panels=2){
  print(Plotname)
  grid.newpage()
  png(paste0(pathToSavePlots,Plotname) ,width=18.5,height=18.5*vertical_panels/2,units="cm",res=600)
  grid.draw(plot)
  dev.off()
}

fmt_dcimals <- function(decimals=0){
  function(x) format(x,nsmall = decimals,scientific = FALSE)
}


# Producing a Sankey flow diagram from the longitudinal data
flowdiagram <- function(initial,final,categories,Age=NULL,reqpal,plottitle,legendtitle){
  ncat <- length(categories)
  if(is.na(which(final == 'NA-In use' | final=='NA-Not in use')[1])){
    # Categories of use only 
    initial <- factor(initial,levels=categories)
    final <- factor(final,levels=categories)
    # Calculate the transition matrix, survey round dependent if required    
    get_tr <- function(tstep){
      if(is.null(Age)){
        trTable <- table(initial,final)
      } else {
        trTable <- table(initial[Age==6*tstep],final[Age==6*tstep])
      }
      tr <- as.matrix(prop.table(trTable,margin=1))
      tr <- ifelse(is.na(tr),0,tr)
      tr[ncat,] <- c(rep(0,ncat-1),1)
    return(tr)}
  } else {
    # Categories of use cross-classified with damage
    extended_cats <- c(categories,'NA-In use','NA-Not in use')
    initial <- factor(initial,levels=extended_cats)
    final <- factor(final,levels=extended_cats)
    #final <- as.factor(ifelse(final == 'NA-In use' & 
    #                            (initial == 'Damaged-In use' | initial == 'Damaged-Not in use'),'Damaged-In use',final))
    get_tr <- function(tstep){
      # Use the overall transition matrix if the Markov assumption holds 
      if(is.null(Age)){
#        tab <- trTable <- table(initial,final)
         tab <- trTable <- tabulate_corrected_categories(initial=initial,final=final,repairs_recoded='none', reweighting=reweighting)
      } else {
         tab <- trTable <- tabulate_corrected_categories(initial=initial[Age==6*tstep],final=final[Age==6*tstep],repairs_recoded='none', reweighting=reweighting)
      }
      trTable[,'Damaged-In use']<- tab[,'Damaged-In use'] * (1 + tab[,'NA-In use']/(tab[,'Damaged-In use']+tab[,'Undamaged-In use'])) 
      trTable[,'Undamaged-In use']<- tab[,'Undamaged-In use'] * (1 + tab[,'NA-In use']/(tab[,'Damaged-In use']+tab[,'Undamaged-In use'])) 
      trTable[,'Damaged-Not in use']<- tab[,'Damaged-Not in use'] * (1 + tab[,'NA-Not in use']/(tab[,'Damaged-Not in use']+tab[,'Undamaged-Not in use'])) 
      trTable[,'Undamaged-Not in use']<- tab[,'Undamaged-Not in use'] * (1 + tab[,'NA-Not in use']/(tab[,'Damaged-Not in use']+tab[,'Undamaged-Not in use'])) 
      trTable1 <- trTable[order(hole_categories),]
      dim <- length(extended_cats)
      tr0 <- matrix(rep(0,dim*dim),nrow=dim,dimnames=list(extended_cats,extended_cats))
      trTable1 <- ifelse(is.na(trTable),tr0,trTable)[1:(dim-2),1:(dim-2)]
      tr <- as.matrix(prop.table(trTable1,margin=1))
      tr <- ifelse(is.na(tr),0,tr)
      tr[ncat,] <- c(rep(0,ncat-1),1)
    return(tr)
  }}
  
  p <- matrix(nrow=9,ncol=ncat)
  
  # Calculate averages for Markov Process
  p[1,] <- c(1,rep(0,ncat-1))
  for(tstep in 1:8){
    tr <- get_tr(tstep)
    p[tstep+1,] <- t(tr) %*% p[tstep,]
  }
  dfp <- as.data.frame(p)
  names(dfp) = categories
  status <- freq <- month <- id <- c()
  idval <- 0
  for(tstep in 1:8){
    tr <- get_tr(tstep)
    for(s1 in 1:ncat){
      for(s2 in 1:ncat){
        idval <- idval+1
        id <- c(id,idval,idval)
        month <- c(month,6*(tstep-1))
        # s1 is initial status (columns in tr matrix)
        status <- c(status,names(dfp)[s1])
        f <- p[tstep,s1]*tr[s1,s2]
        freq <- c(freq,f)
        month <- c(month,6*(tstep)-0.1)
        status <- c(status,names(dfp)[s2])
        freq <- c(freq,f)
        # add duplicate records at the bottom of the file
        if(tstep == 8){
          id <- c(id,idval)
          month <- c(month,6*(tstep))
          status <- c(status,names(dfp)[s2])
          freq <- c(freq,f)
        }
      }
    }
  }
  toPlot <- data.frame(id,month=month,status=as.factor(status),freq)  
  toPlot$status = factor(toPlot$status,levels=categories) 
  library(viridis)
  library(ggplot2)  
  library(ggalluvial)
  p1 <- ggplot(toPlot[toPlot$freq!=0,],
         aes(x = month, stratum = status, alluvium = id,
             y = freq,
             fill = status, label = freq)) +
    geom_alluvium(alpha=0.7)+
    theme_minimal()+
    scale_fill_manual(values=reqpal,name = legendtitle) +
    geom_stratum(alpha=0.8)+
    scale_x_continuous(name = 'Months since distribution',limits=c(0,48.1),breaks = c(0,6,12,18,24,30,36,42,48)) +
    scale_y_continuous(name = 'Proportion of nets')+
    theme_bw()
   savePlot(p1,plottitle,vertical_panels=1)
}


tabulate_corrected_categories = function(initial,final,repairs_recoded, reweighting =reweighting){
  # Summarise input data to provide values for use in model fitting version 
  
  get_weights_for_hole_transitions <- function(initial=initial,final=final){
    # reweight the data matrix for hole categories so that the proportions in the use categories are retained
    
    collapse_use = function(useVar){
      collapsedVar =  factor(recode_factor(useVar,'New'='New','Used last night'='In use','Not used last night'='Not in use','Not in use'='Not in use',
                                    'Damaged-In use'='In use','Damaged-Not in use'='Not in use','Undamaged-In use'='In use','Undamaged-Not in use'='Not in use',
                                    'NA-In use'='In use','NA-Not in use'='Not in use','Attrited'='Attrited'),levels=c("New","In use","Not in use","Attrited"))
    return(collapsedVar)}

    initialu <- collapse_use(transitions$Status.i)
    finalu <- collapse_use(transitions$Status.f)
    hole_initialu <- collapse_use(initial)
    hole_finalu <- collapse_use(final)
    totals_all <- table(initialu,finalu)
    totals_holes <- table(hole_initialu,hole_finalu)
    reweightTot <- totals_all/totals_holes
    reweights <- reweightTot[cbind(as.numeric(hole_initialu),as.numeric(hole_finalu))]
    reweights <- reweights * length(reweights)/sum(reweights)
  return(reweights)}
  
  wts = get_weights_for_hole_transitions(initial=initial,final=final)
  
  if (repairs_recoded != 'none'){
    initial<- recode_factor(initial,'New'='New','Undamaged-Not in use'='S1',
                            'Damaged-Not in use'='S2','Undamaged-In use'='S3','Damaged-In use'='S4')
    final=recode_factor(final,'Undamaged-Not in use'='S1',
                        'Damaged-Not in use'='S2','Undamaged-In use'='S3','Damaged-In use'='S4',
                        'Attrited'='S5','NA-In use'='S6','NA-Not in use'='S7')
    
    # if treating  inconsistent holes as present
    if (repairs_recoded == 'intact'){
      
      # Recode final status to holed, where nets appear to have been repaired
      final[initial=='S4' & final=='S3'] <- 'S4'
      final[initial=='S4' & final=='S1'] <- 'S2'
      final[initial=='S2' & final=='S3'] <- 'S4'
      final[initial=='S2' & final=='S1'] <- 'S2'
      # Recode final status to holed if initial status is holed and final integrity not assessed
      final[initial=='S4' & final=='S6']<-'S4'
      final[initial=='S4' & final=='S7']<-'S2'
      final[initial=='S2' & final=='S6']<-'S4'
      final[initial=='S2' & final=='S7']<-'S2'
    }  
    
    # if treating inconsistent holes as absent
    if (repairs_recoded == 'holed'){
      # Recode initial status to intact, where nets appear to have been repaired
      initial[initial=='S4' & final=='S3'] <- 'S3'
      initial[initial=='S4' & final=='S1'] <- 'S1'
      initial[initial=='S2' & final=='S3'] <- 'S3'
      initial[initial=='S2' & final=='S1'] <- 'S1'
      # Recode final status to holed if initial status is holed and final integrity not assessed
      
      final[initial=='S4' & final=='S6']<-'S4'
      final[initial=='S4' & final=='S7']<-'S2'
      final[initial=='S2' & final=='S6']<-'S4'
      final[initial=='S2' & final=='S7']<-'S2'
    } 
  }
  if(reweighting=='UseWts'){
    data_table <- round(xtabs(wts ~ initial + final)) 
  } else { 
    data_table <- table(initial,final)}
return (data_table)}

select_records <- function(df){
  if(required_net_type != 'All'){
    df <- df[df$Brand==required_net_type,]
  }
  if(required_country != 'All'){
    df <- df[df$Country==required_country,]
  }
  return(df)}

####################### ANALYSIS OF ODE SOLUTIONS #######################
# Forward simulation with Forward Euler and plot 5 year status (including counterfactuals) 
forwardsimulate = function(simoptions=codeMissing,Parameters=testParameters,years=20,tsteps_per_year=1000,
                           plottitle='test.png',requirePlot=TRUE,recycling=FALSE){
  getDifferentials<-function(Parameters){
    qmatrix <- with(Parameters, matrix(data = c(-(u_1+h_1+a_1), 0,  v_3, 0,0,
                                                h_1,-(u_2+a_2),0, v_4,0,  
                                                u_1,0,-(v_3+h_3+a_3),0,0,
                                                0, u_2, h_3, -(v_4+a_4),0,
                                                a_1, a_2, a_3, a_4,0), nrow = 5, ncol = 5, byrow = FALSE,
                                       dimnames = list(c('S_1','S_2','S_3','S_4','A=0'),c('S_1','S_2','S_3','S_4','A')))) 
    return(qmatrix)}
  
  # If recycling is modelled then recycled nets are added into categories S_1 and S_2
    recycle <- function(pstatus=pstatus,qmatrix=qmatrix,status=status,tstep=tstep){
      r <- recyclingModel()
      attrited <- as.numeric(status[5,1]-pstatus[5])
      attrited_with_holes <- as.numeric((t(qmatrix*tstep) %*% (pstatus*c(0,1,0,1,0)))[5])
      # Assume a constant odds ratio 
      
      H <- as.numeric(attrited_with_holes/attrited)
      beta = -r$psi*H - r$psi*r$P3 - 1 + H + r$P3
      alpha = r$psi - 1
      gam= H*r$psi*r$P3
      overall_recycled = r$P3 * attrited
      x = (-beta - sqrt(beta*beta - 4*alpha*gam))/(2*alpha)
      # probability holed net is recycled =  x/H
      recycled_with_holes = ifelse(H > 0, attrited_with_holes * x/H, 0)
      recycled_intact = overall_recycled - recycled_with_holes
      status[1,1] <- status[1,1] + recycled_intact
      status[2,1] <- status[2,1] + recycled_with_holes
      status[5,1] <- status[5,1] - recycled_intact - recycled_with_holes
      recycle_list <- list(status=status, recycled_intact=recycled_intact, recycled_with_holes=recycled_with_holes)
    return(recycle_list)
    }
  
    # used for checking the solution of the quadratic
    checkx <- function(){
      cx = H - x
      bx = r$P3 - x
      dx = 1 - H - r$P3 + x
      sum = x + bx + cx + dx
      psi_check <- x*dx/(bx*cx)  
   return(list(sum,psi_check))}
    
  get_forwardEuler <- function(simoptions,qmatrix,istatus,years,tsteps_per_year,recycling){
    tstep <- 1/tsteps_per_year
    status <- matrix(NA,nrow=5,ncol=1)  
    status = c(as.numeric(istatus))
    statMatrix = matrix(NA,nrow=(years*tsteps_per_year+1),ncol=5)
    statMatrix[1,]=status
    time=1
    sum_attrited_with_holes = 0
    for (year in 1:years){
      for (st in 1:tsteps_per_year){
        time=time+1
        pstatus = status
        status = status + (t(qmatrix*tstep) %*% status)
        sum_attrited_with_holes <- sum_attrited_with_holes + as.numeric((t(qmatrix*tstep) %*% (pstatus*c(0,1,0,1,0)))[5])
        if(recycling) {
          recycle_list = recycle(pstatus=pstatus,qmatrix=qmatrix,status=status,tstep=tstep)
          status = recycle_list$status
          sum_attrited_with_holes <- sum_attrited_with_holes - recycle_list$recycled_with_holes
        }
        statMatrix[time,]=status   
      }
    }
    df <- as.data.frame(statMatrix)
    print(paste0(simoptions, 'attrited_holed: ',sum_attrited_with_holes,' attrited: ',df[nrow(df),5]))   
    colnames(df)=names(istatus)
    df$Month <- as.numeric(seq(1:nrow(df))/tsteps_per_year*12)
    return(df)
  }
  istatus <- with(Parameters,list(S_1=1-U_N, S_2=0, S_3=U_N, S_4=0, A = 0))
  # Uses input parameters
  qmatrix <- getDifferentials(Parameters)  
  estimates <- get_forwardEuler(simoptions=simoptions,qmatrix=qmatrix,istatus,years,tsteps_per_year,recycling)
  # counterfactual1 is no hole acquisition
  # To test the effect of physical robustness uses counterfactual with h_1 =0, h_3 =0
  simoptions2=paste0(simoptions,'_B')
  c1Parameters <- Parameters
  c1Parameters$h_1 = c1Parameters$h_3 =0
  c1matrix <- getDifferentials(c1Parameters)
  counterfactual1 <- get_forwardEuler(simoptions=simoptions2,qmatrix=c1matrix,istatus=istatus,years,tsteps_per_year,recycling)
  # To test the effect of use on holes, specify full use, and prevent transistions to non-use:
  # no use counterfactual is commented out
  simoptions3=paste0(simoptions,'_C')
  c2Parameters <- Parameters
  c2Parameters$U_N = 1 
  #c2Parameters$u_1 = c2Parameters$u_2 =0
  c2Parameters$v_3 = c2Parameters$v_4 =0
  c2matrix <- getDifferentials(c2Parameters)
  #istatus2 <- with(Parameters,list(S_1=1, S_2=0, S_3=0, S_4=0, A = 0))
  istatus2 <- with(Parameters,list(S_1=0, S_2=0, S_3=1, S_4=0, A = 0))
  counterfactual2 <- get_forwardEuler(simoptions=simoptions3,qmatrix=c2matrix,istatus=istatus2,years,tsteps_per_year,recycling)
  estimates$model='A'
  counterfactual1$model='B'
  counterfactual2$model='C'
  allmodels <- rbind(estimates,counterfactual1,counterfactual2)
  # Reorder to match the order of categories in the plot. Better to renumber the parameters.
  allmodels <- allmodels[,c("S_1","S_3","S_4","S_2","A","model","Month")]
  if(requirePlot){
#    facet_bounds <- data.frame(variable=variables,ymin=c(0,0,0,0,0,0),ymax=c(1,1,1,1,1,1),x=rep(6,6))
    pdata <- melt(allmodels,id=c("model","Month"))
    names(pdata) = c("model",'Month','Category','Prop')
    levels(pdata$Category)= hole_categories[2:6]
    reqpal=modifiedinferno[2:6]
    p1 = ggplot(pdata,
                aes(x = Month, stratum = Category, y = Prop,
                    fill = Category)) +
      geom_area(stat="identity", alpha=0.7)+
      theme_minimal()+
      scale_fill_manual(values=reqpal) +
      scale_x_continuous(name = 'Months since distribution',limits=c(0,48),breaks = c(0,6,12,18,24,30,36,42,48)) +
      scale_y_continuous(name = 'Proportion of nets')+
      facet_wrap(~ model)+
      theme_bw()
    savePlot(p1,plottitle,vertical_panels=0.6)
  }
  return(allmodels)}


get_propAttrition_in_holed <- function(e=allmodels,p=Parameters,tsteps_per_year=tsteps_per_year){
    e <- e[e$model=='A',]
    attrition <- e$A - c(0,e$A[1:(nrow(e)-1)]) 
    prop_in_holed <- (p$a_2*e$S_2 +  p$a_4*e$S_4)/(p$a_1*e$S_1 +  p$a_2*e$S_2 + p$a_3*e$S_3 +  p$a_4*e$S_4)
    prop_holed <- e$S_2 + e$S_4
    attrition_in_holed <- attrition*prop_in_holed
    lifetime_holed <- sum(prop_holed)/tsteps_per_year
    return(list(pAttrition_in_holed=sum(attrition_in_holed)/sum(attrition),lifetime_holed=lifetime_holed))
}

get_interval_estimates <- function(samples_posterior=samples_posterior,plottitle=plottitle){
  useLastNight <- get_use_last_night()
  pars <- samples_posterior[sample(nrow(samples_posterior), nsamples_posterior), ] 
  for (i in 1:nsamples_posterior){
    print(i)
    h_1 <- pars$h1[i]
    h_3 <- pars$h3[i]
    u_1 <- pars$u1[i]
    u_2 <- pars$u2[i]
    v_3 <- pars$v3[i]
    v_4 <- pars$v4[i]
    a_1 <- pars$a1[i]
    a_2 <- pars$a2[i]
    a_3 <- pars$a3[i]
    a_4 <- pars$a4[i]
    P_0 <- pars$p0[i]
    U_N <- 0
    U_H <- useLastNight$U_H
    U_I <- useLastNight$U_I
    parameters <- list(h_1=h_1 , h_3=h_3, u_1=u_1, u_2=u_2, 
                       v_3=v_3 , v_4=v_4, a_1=a_1 , a_2=a_2, 
                       a_3=a_3 , a_4=a_4, P_0=P_0 , U_N=U_N, U_H=U_H , U_I=U_I)
    stats <- get_summaryStats(simoptions=NULL,i=i, Parameters=parameters,plottitle='null',requirePlot=FALSE)
    if(i==1){
      df <- as.data.frame(list(parameters,stats))
    } else {
      statsdf <- as.data.frame(list(parameters,stats))
      df <- rbind(df,statsdf)
    }
  }
  interval_estimates <- as.data.frame(as.list(c(names(df)[1],quantile(df[,1], probs = c(0.025,0.5,0.975),na.rm = TRUE))))
  names(interval_estimates)[1] <- 'Variable'
  for(var in c(2:14,24:37)){
    df1 <- as.data.frame(as.list(c(names(df)[var],quantile(df[,var], probs = c(0.025,0.5,0.975),na.rm = TRUE))))
    names(df1)[1] <- 'Variable'
    interval_estimates <- rbind(interval_estimates,df1)
  }
  interval_estimates$model <- 'A'
  for (model in c('A','B','C')){
    df2 <- df[which(df$summaryByModel.models==model),]
    for(var in c(16:23)){
      df1 <- as.data.frame(as.list(c(names(df2)[var],quantile(df2[,var], probs = c(0.025,0.5,0.975),na.rm = TRUE))))
      names(df1)[1] <- 'Variable'
      df1$model <- model
      interval_estimates <- rbind(interval_estimates,df1)
    }
  }
return(interval_estimates)}  
  
# calculate summary statistics and plots from estimates of the ODE parameters 
get_summaryStats <- function(simoptions,i, Parameters=Parameters,years=20,tsteps_per_year=1000,plottitle=NULL,requirePlot=FALSE){
  simSets <- forwardsimulate(simoptions=simoptions,Parameters=Parameters,years=20,tsteps_per_year=1000,plottitle=plottitle,requirePlot=requirePlot,recycling=FALSE)
  models <- c('A','B','C')
  # A: observed; B: no holes; C: 100% use
  propLifetimeHoled <- propNotAttrited <- propInUse <- propInUseLastNight <- propHoled <- rep(NA,3)
  medianLife <- meanLife <- propLifetimeInUse <-propLifetimeInUseLastNight <- rep(NA,3)
  
  for (m in 1:3){
    subSet <- simSets[simSets$model==models[m],] 
    #  Average life of LLIN (years)
    meanLife[m] <- sum(1-subSet$A)*subSet$Month[1]/12
    medianLife[m] <- subSet[subSet$A > 0.5,]$Month[1]/12
    # Proportion of simulation period for which net is not attrited  
    propNotAttrited[m] <- sum(1-subSet$A)/length(subSet$A)
    # Proportion of simulation period for which net is in use  
    propInUse[m] <- with(subSet, (sum(S_3) + sum(S_4))/length(A))
    # Proportion of simulation period for which net is in use lastnight 
    propInUseLastNight[m] <- with(subSet, (sum(S_3*Parameters$U_I) + sum(S_4*Parameters$U_H))/length(A))  
    # Proportion of simulation period for which net is holed  
    propHoled[m] <- with(subSet, (sum(S_2) + sum(S_4))/length(A))
    # Proportion of lifetime of LLIN for which it is in 'use'
    propLifetimeInUse[m] <- with(subSet,(sum(S_3) + sum(S_4))/sum(1-A))
    # Proportion of lifetime of LLIN for which it is holed
    propLifetimeHoled[m] <- with(subSet,(sum(S_2) + sum(S_4))/sum(1-A))  
    # Proportion of lifetime of LLIN for which was in use last night i.e. indeed used  
    propLifetimeInUseLastNight[m] <- with(subSet,(sum(S_3*Parameters$U_I) + sum(S_4*Parameters$U_H))/sum(1-A))
  }
  print(paste('Mean lifetimes',meanLife[1],meanLife[2],meanLife[3]))
  # Reduction in net lifetime attributable to holes (years) 
  reductionInLifetimeDueToHoles <- medianLife[2] - medianLife[1]
  propLossLifetimeDueToHoles <- 1 - propNotAttrited[1]/propNotAttrited[2] 
  # Proportion of lack of use attributable to holes (comparison with counterfactual)
  propLackOfUseDueToHoles <- 1 - propInUse[1]/propInUse[2]
  propLackOfUseLastNightDueToHoles <- 1 - propInUseLastNight[1]/propInUseLastNight[2]
  RelativeRateHolesByUse <- Parameters$h_3/Parameters$h_1
  #Proportion of attrition in holed nets and lifetime of holed nets
  propAttrition_in_holed <- get_propAttrition_in_holed(e=simSets,p=Parameters,tsteps_per_year=tsteps_per_year)
  propAttritionInHoled <- propAttrition_in_holed$pAttrition_in_holed
  medianLifeHoled <- propAttrition_in_holed$lifetime_holed
  totalNightsInUse <- 365*medianLife[1]*propLifetimeInUseLastNight[1]
  totalNightsInUseCounterfactual <- 365*medianLife[2]*propLifetimeInUseLastNight[2]
  propNightsLost <- 1-  totalNightsInUse/totalNightsInUseCounterfactual
  # Use P_V = 0.08 from Briet et al 2020
  P_V =0.08
  propImpactLost <- 1 - (1-P_V )*(1- propNightsLost)
  RelativeRateUseByHoles <- Parameters$u_2/Parameters$u_1

  summaryByModel <- data.frame(models,propLifetimeHoled,propNotAttrited,propInUse,propInUseLastNight, 
                               propHoled,medianLife,propLifetimeInUse,propLifetimeInUseLastNight,
                               propAttritionInHoled,medianLifeHoled)
  summaryStats <- list(summaryByModel=summaryByModel,
                       reductionInLifetimeDueToHoles=reductionInLifetimeDueToHoles,
                       propLossLifetimeDueToHoles=propLossLifetimeDueToHoles,
                       propLackOfUseDueToHoles=propLackOfUseDueToHoles,
                       propLackOfUseLastNightDueToHoles=propLackOfUseLastNightDueToHoles,
                       RelativeRateHolesByUse=RelativeRateHolesByUse,
                       propAttritionInHoled=propAttritionInHoled,
                       medianLifeHoled=medianLifeHoled,
                       totalNightsInUse=totalNightsInUse,
                       totalNightsInUseCounterfactual=totalNightsInUseCounterfactual,
                       propNightsLost=propNightsLost,
                       propImpactLost=propImpactLost,
                       RelativeRateUseByHoles=RelativeRateUseByHoles)
return(summaryStats)}

# calculate use last night in the overall dataset as a proportion of all used nets, by hole status
get_use_last_night <- function(){
  df <- read_pooled_llin_data()
  tab <- with(df[!is.na(df$PHI),],table(Status,ifelse(PHI>0,'holed','intact')))
  U_H <- tab[2,1]/(tab[2,1]+tab[3,1])
  U_I <- tab[2,2]/(tab[2,2]+tab[3,2])
  uln <- list(U_H=U_H,U_I=U_I)
return(uln)}

# extract parameters from data frame
get_parameters = function(samples_posterior){
  pln <- get_use_last_night()
  h_1 <- median(samples_posterior$h1)
  h_3 <- median(samples_posterior$h3)
  u_1 <- median(samples_posterior$u1)
  u_2 <- median(samples_posterior$u2)
  v_3 <- median(samples_posterior$v3)
  v_4 <- median(samples_posterior$v4)
  a_1 <- median(samples_posterior$a1)
  a_2 <- median(samples_posterior$a2)
  a_3 <- median(samples_posterior$a3)
  a_4 <- median(samples_posterior$a4)
  P_0 <- median(samples_posterior$p0)
  U_N <- 0
  U_H <- pln$U_H
  U_I <- pln$U_I
  parameters <- list(h_1=h_1 , h_3=h_3, u_1=u_1, u_2=u_2, 
                     v_3=v_3 , v_4=v_4, a_1=a_1 , a_2=a_2, 
                     a_3=a_3 , a_4=a_4, P_0=P_0 , U_N=U_N, U_H=U_H , U_I=U_I)
  return(parameters)}

# Parameters are
# h_1	Acquisition of holes in unused nets
# h_3	Acquisition of holes in nets in use
# u_1	Putting intact nets to use
# u_2	Putting holed nets to use
# v_3	Taking intact nets out of use
# v_4	Taking holed nets out of use
# a_1	Attrition of intact, unused nets
# a_2	Attrition of holed, unused nets
# a_3	Attrition of intact, used nets
# a_4	Attrition of holed, used nets
# P_0	Proportion of nets at follow-up for which physical integrity was evaluated 
# U_N	Probability that a new net is taken into use immediately on receipt
# U_H	Probability that a holed net was used on night before survey
# U_I	Probability that an intact net was used on night before survey

# Overall workflow for ODE models.
processODE <- function(required_net_type,required_country){

  # if the interval estimates are not estimated then the function returns NULL
  interval_estimates = NULL
  if (requireDataDescription) {
    if(nrow(dfw)>0 & required_net_type=='All'){
    # get table for analysis of whether initial hole index predicts transitions
    PHItable <- get_PHItable(transitions=transitions)
    
    # Flow diagram for observed transition probabilities time independent probabilities  
    flowdiagram(initial=transitions$Status.i,final=transitions$Status.f,
                categories=categories,reqpal = modified,
                plottitle=paste0('FlowConstantUseCategories',required_options(),'.png'),
                legendtitle='Use categories')
    # Flow diagram for observed transition probabilities time dependent probabilities  
    flowdiagram(initial=transitions$Status.i,final=transitions$Status.f,categories=categories,Age=transitions$Age.f,
                reqpal = modified,
                plottitle=paste0('FlowAgeDependentUseCategories',required_options(),'.png'),
                legendtitle='Use categories')
    # Flow diagram for observed transition probabilities (including repairs) assuming Markov process 
    flowdiagram(initial=hole_transitions$HoleStatus.i,final=hole_transitions$HoleStatus.f,
                categories=hole_categories,reqpal=modifiedinferno,
                plottitle=paste0('FlowConstantHoleCategories',required_options(),'.png'),
                legendtitle='Use/damage categories')
    # Flow diagram for observed transition probabilities time dependent probabilities 
    flowdiagram(initial=hole_transitions$HoleStatus.i,final=hole_transitions$HoleStatus.f,
                categories=hole_categories,Age=hole_transitions$Age.f,reqpal = modifiedinferno,
                plottitle=paste0('FlowAgeDependentHoleCategories',required_options(),'.png'),
                legendtitle='Use/damage categories')
    # Save dataset for fitting of ODE model
    filename <- paste0('hole_transitions',required_options(),'.csv')
    write.csv(hole_transitions,file=filename)
    } 
  }  
  if (requireODEFitting){
    if (nrow(dfw)>0) {
      source('../MarkovmodelStan/netdurability_fixqp.R')
      if (repairs_recoded == 'intact'){
        data_matrix <- matrix(tabulate_corrected_categories(initial=hole_transitions$HoleStatus.i,final=hole_transitions$HoleStatus.f,
          repairs_recoded = 'intact', reweighting =reweighting),nrow=5,ncol=7)
        print(paste0('*i: ',required_country,' ',required_net_type))
        samples_posterior <- netdurability_run_fixqp(counts=data_matrix)
        setwd('./PMI_Markov')
        write.csv(samples_posterior,file=paste0('samples_posterior',required_options(),'.csv'))
      }
      if (repairs_recoded == 'holed'){
        data_matrix <- matrix(tabulate_corrected_categories(initial=hole_transitions$HoleStatus.i,final=hole_transitions$HoleStatus.f,
          repairs_recoded = 'holed', reweighting =reweighting),nrow=5,ncol=7)
        print(paste0('*h: ',required_country,' ',required_net_type))
        samples_posterior <- netdurability_run_fixqp(counts=data_matrix)
        setwd('./PMI_Markov')
        write.csv(samples_posterior,file=paste0('samples_posterior',required_options(),'.csv'))
      }
    }
  }
  if(requireODEpostprocessing){
    plottitle <-   paste0('Predictions_',required_options(),'.png')
    #Use as reference the analysis treating repaired nets as holed
    filename=paste0('samples_posterior',required_options(),'.csv')
    if(file.exists(filename)){
      samples_posterior <- read.csv(filename)
      pars <- get_parameters(samples_posterior)
      summaryStats <- get_summaryStats(simoptions=NULL,i=NULL,Parameters=pars,years=10,tsteps_per_year=1000,plottitle=plottitle,requirePlot=TRUE)
    }
    if(requireCalculationIntervalEstimates){
      interval_estimates <- get_interval_estimates(samples_posterior=samples_posterior,plottitle=plottitle)
      interval_estimates$net_type <- required_net_type
      interval_estimates$country <- required_country
      interval_estimates$repair_recoding = repairs_recoded
    }
  }
  if(!is.null(interval_estimates)){
    write.csv(interval_estimates,file=paste0('interval_estimates',required_options(),'.csv'))
  } else { 
    interval_estimates <- read.csv(file=paste0('interval_estimates',required_options(),'.csv'))
    interval_estimates$X = NULL
  }
  interval_estimates$X50. <- as.numeric(as.character(interval_estimates$X50.))
  interval_estimates$X2.5. <- as.numeric(as.character(interval_estimates$X2.5.))
  interval_estimates$X97.5. <- as.numeric(as.character(interval_estimates$X97.5.))
  return(interval_estimates)
}

# To plot country comparisons of ODE outputs  
plotComparisons <- function(all_interval_estimates,classifier,filterval){
  if(classifier=='net_type'){
    filtervar <- all_interval_estimates$country
    plottitle = paste0('Comparison_of_net_types_',filterval,codeMissing,repairs_recoded,reweighting,'.png')
    all_interval_estimates$category <- all_interval_estimates$net_type
    legendTitle = 'Net type'
  } else { 
    filtervar <- all_interval_estimates$net_type
    plottitle = paste0('Comparison_of_countries_',filterval,codeMissing,repairs_recoded,reweighting,'.png')
    all_interval_estimates$category <- all_interval_estimates$country
    legendTitle = 'Country'
  }

  plotDataL <- with(all_interval_estimates,all_interval_estimates[filtervar==filterval & model == 'A' & repair_recoding == repairs_recoded &
                                                                    !is.na(match(Variable,c("summaryByModel.medianLife","reductionInLifetimeDueToHoles"))),])
  plotDataM <- with(all_interval_estimates,all_interval_estimates[filtervar==filterval & model == 'A' & repair_recoding == repairs_recoded &
                                                                    !is.na(match(Variable,c("h_1","h_3","u_1","u_2","v_3","v_4","a_1","a_2","a_3","a_4"))),])
  plotDataP <- with(all_interval_estimates,all_interval_estimates[filtervar==filterval & model == 'A' & repair_recoding == repairs_recoded &
                                                                    !is.na(match(Variable,c("P_0",
                                                                                            "summaryByModel.propLifetimeInUse",
                                                                                            "summaryByModel.propLifetimeInUseLastNight",
                                                                                            "propLackOfUseDueToHoles",
                                                                                            "propLackOfUseLastNightDueToHoles",
                                                                                            "propLossLifetimeDueToHoles"))),])
  plotDataL$Variable <- factor(plotDataL$Variable,
                               levels = c("summaryByModel.medianLife","reductionInLifetimeDueToHoles"))
  plotDataM$Variable <- factor(plotDataM$Variable,
                               levels = c("h_1","h_3","u_1","u_2","v_3","v_4","a_1","a_2","a_3","a_4"))
  plotDataP$Variable <- factor(plotDataP$Variable,
                               levels = c(
                                 "P_0",
                                 "summaryByModel.propLifetimeInUse",
                                 "summaryByModel.propLifetimeInUseLastNight",
                                 "propLackOfUseDueToHoles",
                                 "propLackOfUseLastNightDueToHoles",
                                 "propLossLifetimeDueToHoles"))
 
  
  
  library(cowplot)
  #size - (default: 0.5) thickness of the lines
  #linetype - (default: 1=solid) the type of the lines
  #colour - (default: "black") the color of the lines
  #width - (default: 0.9) width of the whiskers
  #alpha - (default: 1=opaque) the transparency of the lines 
  p3 <- ggplot(parse = TRUE)  + theme_bw() + scale_y_continuous(name = 'Duration (years)')+ 
    geom_errorbar(data=plotDataL[plotDataL$category=='All',], aes(x = Variable, ymax=as.numeric(X97.5.),ymin=as.numeric(X2.5.)),size=2,width=0,colour='grey',alpha=0.8) + 
    geom_point(data=plotDataL, aes(x = Variable, colour = category,
                                   y = as.numeric(X50.)), show.legend=TRUE) + 
    scale_colour_discrete(name=legendTitle) +
    geom_point(data=plotDataL[plotDataL$category=='All',], aes(x = Variable, y = as.numeric(X50.)),colour='black',size=2, shape=15) + 
    scale_x_discrete(name=element_blank(), labels= c("Median Lifetime (years)","Reduction in Lifetime"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)
    )
  p1 <- ggplot(parse = TRUE) + theme_bw() + scale_y_log10(name = 'Rate (per year)') +
    geom_errorbar(data=plotDataM[plotDataM$category=='All',], aes(x = Variable, ymax=as.numeric(X97.5.),ymin=as.numeric(X2.5.)),size=2,width=0,colour='grey',alpha=0.8) + 
    geom_point(data=plotDataM, aes(x = Variable, colour = category,
                                   y = as.numeric(X50.)), show.legend=FALSE) + 
    geom_point(data=plotDataM[plotDataM$category=='All',], aes(x = Variable, y = as.numeric(X50.)),colour='black',size=2, shape=15) + 
    scale_x_discrete(name=element_blank(),labels= c(expression(h[1]),expression(h[3]),
                                                    expression(u[1]),expression(u[2]),expression(v[3]),expression(v[4]),expression(a[1]),
                                                    expression(a[2]),expression(a[3]),expression(a[4])))
  p2 <- ggplot(parse = TRUE)  + theme_bw() + scale_y_continuous(name = 'Proportion',limits = c(0, 1)) +
    geom_errorbar(data=plotDataP[plotDataP$category=='All',], aes(x = Variable, ymax=as.numeric(X97.5.),ymin=as.numeric(X2.5.)),size=2,width=0,colour='grey',alpha=0.8) + 
    geom_point(data=plotDataP, aes(x = Variable, colour = category,
                                   y = as.numeric(X50.)), show.legend=FALSE) + 
    geom_point(data=plotDataP[plotDataP$category=='All',], aes(x = Variable, y = as.numeric(X50.)),colour='black',size=2, shape=15) + 
    scale_x_discrete(name=element_blank(),labels= c(expression(paste(P[0] )),
                                                    expression("Lifetime in Use"^{"1"}),expression("Lifetime in Use"^{"2"}),
                                                    expression("Non-Use due to holes"^{"1"}),expression("Non-Use due to holes"^{"2"}),"Loss in Lifetime due to holes")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) 
  p4 <- plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12,ncol=3)
  savePlot(p4,Plotname=plottitle,vertical_panels=1)
}

plotPHIscatter = function(required_country='Kenya'){
  dfw <- select_records(DataWide)
  transitions <- get_transitions(dfw=dfw)
  hole_transitions <- get_hole_transitions(transitions)
  ppHI = ggplot(data=hole_transitions, aes(x=log(PHI.i+1),y=log(PHI.f+1)))+  theme_bw()+
    theme(text = element_text(size = 20)) +
    geom_jitter(position = position_jitter(width = 0.2, height = 0.2),shape=1) +
    geom_abline(intercept = 0, slope = 1) +
    geom_vline(xintercept = log(21),linetype='dashed') +
    geom_hline(yintercept = log(21),linetype='dashed') +
    scale_x_continuous(name='Log(initial pHI +1)')+
    scale_y_continuous(name='Log(final pHI +1)')
  savePlot(ppHI,Plotname=paste0('PHIScatter.png'),vertical_panels=2)
}



recyclingModel <- function(P1 = 0.927891982,
                           P2 = 0.631508196928867,
                           P3 = 0.76469,
                           P5 = 0.461049247182706){

# From DXX model:
# P1 Proportion of destroyed nets that were holed before destruction

P1 = 0.771008356683577/0.782814658444883
  
# From DRX model:
# Proportion of destroyed and removed nets that were holed before destruction

P2 = 0.90187444883464/0.997120476949762
  
# From questionnaire responses, proportion removed among removed and destroyed nets

P3 = 5687/(1750+5687) 

# Define P4 as proportion of removed nets that were holed before destruction:
# P2 = P1(1-P3) + P4P3

P4 = (P2 - P1*(1-P3))/P3

# From DRA model:
# Proportion of 'absent' nets that were holed before absence


P5 = 0.834087789924736/0.999919629104698

# Define P6 as the proportion of 'lost' nets that are destroyed, so that (1-P6) is the proportion that survive:
# P5 = P4(1-P6) + P1P6
# Rearranging, this gives:

P6 = (P5 - P4)/(P1 - P4)

# This gives a negative, use 1- P3

P6 = 1 - P3

# Odds ratio used for calculating time step specific proportion holed among recycled
psi = P4*(1-P1)/(P1*(1-P4))

returnList = list(psi=psi,P1=P1,P2=P2,P3=P3,P4=P4,P5=P5,P6=P6)

return(returnList)}

RecyclingAnalysis <- function(){
  filename=paste0('samples_posteriorAllKenyaDXAholedUseWts.csv') 
  # Recycling is plotted with parameter set DXA
  samples_posterior <- read.csv(filename)
  pars <- get_parameters(samples_posterior)
  modelRRR <- forwardsimulate(simoptions='RRR',Parameters=pars,years=30,tsteps_per_year=1000,plottitle='Predictions_AllkenyaRRRholedUseWts.png',requirePlot=TRUE,recycling=TRUE)
  modelRRR <- modelRRR[modelRRR$model=='A',]
  modelRRR$model <- 'A1R'
  listA1 = get_simulationResults(filename='samples_posteriorAllKenyaDXAholedUseWts.csv',modelTitle='A1')
  listA2 = get_simulationResults(filename='samples_posteriorAllKenyaDXXholedUseWts.csv',modelTitle='A2')
  listA3 = get_simulationResults(filename='samples_posteriorAllKenyaDRXholedUseWts.csv',modelTitle='A3')
  listA4 = get_simulationResults(filename='samples_posteriorAllKenyaDRAholedUseWts.csv',modelTitle='A4')
  concatenated_estimates=rbind(listA1$interval_estimates,listA2$interval_estimates,listA3$interval_estimates,listA4$interval_estimates)
  tabulate_interval_estimates(cumulated_interval_estimates=concatenated_estimates,outfile='EstimatesFromRecyclingAnalysis.csv')

  fourModelsToPlot <- rbind(listA2$modeldf,listA3$modeldf,listA4$modeldf,modelRRR)
  pdata <- melt(fourModelsToPlot,id=c("model","Month"))
  plottitle <- 'ForwardSimulationWithRecyclingKenya.png'
  names(pdata) = c("model",'Month','Category','Prop')
  levels(pdata$Category)= hole_categories[2:6]
  reqpal=modifiedinferno[2:6]
    p1 = ggplot(pdata,
                aes(x = Month, stratum = Category, y = Prop,
                    fill = Category)) +
      geom_area(stat="identity", alpha=0.7)+
      theme_minimal()+
      scale_fill_manual(values=reqpal) +
      scale_x_continuous(name = 'Months since distribution',limits=c(0,48),breaks = c(0,6,12,18,24,30,36,42,48)) +
      scale_y_continuous(name = 'Proportion of nets')+
      facet_wrap(~ model)+
      theme_bw()
    savePlot(p1,plottitle,vertical_panels=1.2)
 return(fourModelsToPlot)}
  
# Compare estimates from different definitions of damage
SensitivityAnalysis <- function(){
  listA1 = get_simulationResults(filename=paste0('samples_posteriorAllKenyaDXAholedUseWts.csv'),modelTitle='A1') 
  listA1_I = get_simulationResults(filename='samples_posteriorAllKenyaDXAintactUseWts.csv',modelTitle='A1_I') 
  listA1_0 = get_simulationResults(filename='samples_posteriorAllKenyaDXAholedUseWts0.csv',modelTitle='A1_0')
  concatenated_estimates=rbind(listA1$interval_estimates,listA1_I$interval_estimates,listA1_0$interval_estimates)
  tabulate_interval_estimates(cumulated_interval_estimates=concatenated_estimates,outfile='EstimatesFromSensitivityAnalysis.csv')
  
  threeModelsToPlot <- rbind(listA1$modeldf,listA1_I$modeldf,listA1_0$modeldf)
  pdata <- melt(threeModelsToPlot,id=c("model","Month"))
  plottitle <- 'ForwardSimulationSensitivityKenya.png'
  names(pdata) = c("model",'Month','Category','Prop')
  levels(pdata$Category)= hole_categories[2:6]
  reqpal=modifiedinferno[2:6]
  p1 = ggplot(pdata,
              aes(x = Month, stratum = Category, y = Prop,
                  fill = Category)) +
    geom_area(stat="identity", alpha=0.7)+
    theme_minimal()+
    scale_fill_manual(values=reqpal) +
    scale_x_continuous(name = 'Months since distribution',limits=c(0,48),breaks = c(0,6,12,18,24,30,36,42,48)) +
    scale_y_continuous(name = 'Proportion of nets')+
    facet_wrap(~ model)+
    theme_bw()
  savePlot(p1,plottitle,vertical_panels=0.6)
return(threeModelsToPlot)}

# Read MCMC output file and calculate interval estimates
get_simulationResults = function(filename,modelTitle){
  samples_posterior <- read.csv(filename)
  pars <- get_parameters(samples_posterior)
  modeldf <- forwardsimulate(simoptions=modelTitle,Parameters=pars,years=30,tsteps_per_year=1000,requirePlot=FALSE,recycling=FALSE)
  modeldf <- modeldf[modeldf$model=='A',]
  modeldf$model <- modelTitle
  interval_estimates <- get_interval_estimates(samples_posterior=samples_posterior,plottitle=NULL)
  interval_estimates$grouping=modelTitle
  resultsList = list(modeldf=modeldf,interval_estimates=interval_estimates)
  return(resultsList)
}

# Tabulate parameter estimates
tabulate_interval_estimates = function(cumulated_interval_estimates,outfile){
  toTabulate <- cumulated_interval_estimates[cumulated_interval_estimates$model=='A',]
  toTabulate$X50. = as.numeric(as.character(toTabulate$X50.))
  toTabulate$X2.5. = as.numeric(as.character(toTabulate$X2.5.))
  toTabulate$X97.5. = as.numeric(as.character(toTabulate$X97.5.))
  temp1 <- dcast(na.omit(toTabulate), grouping ~ Variable, fun=mean, value.var='X50.')
  temp2 <- dcast(na.omit(toTabulate), grouping ~ Variable, fun=mean, value.var="X2.5.")
  temp3 <- dcast(na.omit(toTabulate), grouping ~ Variable, fun=mean, value.var="X97.5.")
  temp1$percentile=50
  temp2$percentile=2.5
  temp3$percentile=97.5
  allEstimates = rbind(temp1,temp2,temp3)
  write.csv(allEstimates,file=outfile)
}


required_options = function() {return (paste0(required_net_type,required_country,codeMissing,repairs_recoded,reweighting))}
