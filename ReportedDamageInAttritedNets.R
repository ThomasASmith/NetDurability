getReportedDamageByAge = function(){
  df <- survey_data
  df1 = df[df$Status=='Attrited',]
  df2 = as.data.frame(with(df1,table(Round,Damage_related)))
  df3 = df2 %>%
    group_by(Round,Damage_related) %>%
    summarise(r_t = Freq) %>%
    mutate(prop = r_t/sum(r_t),total=sum(r_t))

  myfxn= function(X,total,arg){ res=binom.test(X,total)
     return(res$conf.int[arg])}
  toPlot = data.frame(df3[df3$Damage_related=='Yes',])
  toPlot$lowerL=mapply(X=toPlot$r_t,FUN=myfxn,total=toPlot$total,arg=1)
  toPlot$upperL=mapply(X=toPlot$r_t,FUN=myfxn,total=toPlot$total,arg=2)
  toPlot$Age.f = as.numeric(as.character(toPlot$Round))
  # add additional records to enable plotting of projections over longer net lifetimes
  additional_records =data.frame(Age.f=c(0.08,60,72,84,96))
  additional_records[setdiff(names(toPlot), names(additional_records))] <- NA
  toPlot =rbind(toPlot,additional_records)
  # fit logistic regression curve in age of net
  model <- glm(cbind(r_t,total-r_t) ~ log(Age.f), data = toPlot, family = binomial)
  lp = model$coefficients[1] + log(toPlot$Age.f)*model$coefficients[2]
  toPlot$predict = 1/(1 + exp(-lp))
  plot = ggplot(data= toPlot, aes(x=Age.f, y = prop))+
    geom_point(stat = "identity",size=5) + geom_errorbar(aes(ymin = lowerL, ymax = upperL),width = 2, position=position_dodge(.9))+
    geom_smooth(aes(y=predict),se=F)+
    scale_x_continuous(name = 'Months since distribution',limits=c(0,100),breaks = c(0,12,24,36,48,60,72,84,96)) +
    scale_y_continuous(name = 'Proportion of attrited nets', limits=c(0,1))+
    theme_bw(base_size = 24)
return(list(analysis=toPlot,plot=plot))}

getReportedDamageByPreviousDamage = function(reported){
  ### Computation of proportion of attrition from damaged nets in projection from ODE
  filename='interval_estimatesAllDXAholedUseWts.csv'
  interval_estimates <- read.csv(filename)
  pars = as.list(interval_estimates$X50.)
  names(pars)=as.character(interval_estimates$Variable)
  projections = read.csv('Projections from model A1.csv')
  projections$p_d=with(projections,(S_2 + S_4)/(1- A))
  projections$attrition_from_S_1=with(projections, S_1*pars$a_1/(1- A))
  projections$a =with(projections,(S_1*pars$a_1 + S_2*pars$a_2 + S_3*pars$a_3 + S_4*pars$a_4)/(1- A))
  projections$halfyear = round((projections$Month+3)/6)
  # average attrition rate and convert to annual rate
  aggregated = aggregate(projections[,9:10],list(projections$halfyear),sum) /500
  aggregated$p_d = aggregate(projections$p_d,list(projections$halfyear),mean)[,2]
  aggregated <- cbind(aggregated[1:8,],reported[1:8,])
  ReportedDamageByPreviousDamage = ggplot(data=aggregated, aes(x=p_d, y = prop))+
    geom_point(stat = "identity",size=5) + geom_errorbar(aes(ymin = lowerL, ymax = upperL),width = 2, position=position_dodge(.9))+
    scale_x_continuous(name = 'Estimated proportion with previous damage',limits=c(0,1)) +
    scale_y_continuous(name = 'Proportion reported as due to damage', limits=c(0,1))+
    geom_smooth(aes(y=predict),se=F)+
    theme_bw(base_size = 14)
  return(list(dataset=aggregated,plot=ReportedDamageByPreviousDamage,a_1=pars$a_1))
}


getModelForTypeOfAttrition = function(input){
    library(rjags)
    library(coda)
    load.module("mix")
    library(runjags)
    load.module("dic")
    load.module("glm")
    dataset=input$dataset
    a_1 = input$a_1

    jagscode = "model{

  a_c0 ~dunif(0,1)
  a_c <- a_c0*a_1

  # the minimum value of a_i corresponds to no incremental damage (a_i = 0)
  # the maximum value of a_i corresponds to a_u = 0
  # hence: a_i = (a -a_c)/p_d at the time point with the minimum ratio a/p_d

  a_i0 ~dunif(0,1)
  a_i <-a_i0*(a[8] - a_c)/p_d[8]

  # Proportion of incrementally damaged nets that attrite with damage reported as cause of attrition
  p_ix ~dunif(0,1)

  for (t in 1:surveys){

    # Proportion of attrited nets catastrophically damaged at point of attrition
    p_c[t] <- a_c/a[t]

    # Proportion of attrited nets incrementally damaged at point of attrition
    p_i[t] <- p_d[t]*a_i/a[t]

    # Proportion of attrited nets undamaged at point of attrition

    p_u[t] <- 1 -  p_c[t] - p_i[t]
    a_u[t] <- p_u[t] * a[t]

    # Proportion of attrited nets in which the owner attributes attribution to damage
    p_x[t] <- (a_c + a_i*p_d[t]*p_ix)/a[t]
    r_x[t] ~ dbin(p_x[t],n_x[t])

    # Proportion of attrited nets previously damaged
    PD0[t] <- p_d[t]*(a_c + a_u[t] + a_i)/a[t]

    # Overall proportion of attrited nets damaged when attrited
    PD1[t] <- p_c[t] + p_i[t]

    C1[t] <- p_x[t]-p_c[t]*(1-p_d[t]) # Previously damaged and attributed to damage
    C2[t] <- p_c[t]*(1-p_d[t])  # Previously undamaged, but damaged at end of life (Always attributed to damage)
    C3[t] <- p_u[t]*p_d[t] + p_i[t] + p_c[t] - p_x[t] # Previously damaged but not attributed to damage
    C4[t] <- p_u[t]*(1-p_d[t])  # Unrelated to damage
    }
  }
  "

  # DATA are:
  datajags <- list(
    p_d = dataset$p_d,
    a = dataset$a,
    a_1=a_1,
    n_x =dataset$total,
    r_x =dataset$r_t,
    # p_d[t]  Proportion of attrition in nets damaged up to point of attrition by time interval
    # a[t]    Total attrition by time interval
    # r_x[t]  Number of attrited nets with damage reported as cause of attrition by survey
    # n_x[t]  Number of attrited nets with data on cause of attrition
    surveys=8)

  inits<-list(list(),list(),list())
  jags.model.out = jags.model(textConnection(jagscode), data=datajags, inits = list(.RNG.name = "base::Wichmann-Hill"), n.adapt=1000, n.chains=3)
  model.samples <- coda.samples(jags.model.out, c("a_c", "a_i","a_u","p_x","p_ix","p_i","p_u","C1","C2","C3","C4"),n.iter = 4000)
  modelResults <- summary(model.samples)
  #plot(model.samples, trace=FALSE, density = TRUE)
  jagsResults=as.data.frame(modelResults$quantiles)
  dataset$p_x = jagsResults[1:8,3]
  dataset$p_x_lower = jagsResults[1:8,1]
  dataset$p_x_upper = jagsResults[1:8,5]
  a_c = jagsResults[33,3]
  cats=c('a: Previous damage, damage reported as reason',
         'b: Previously undamaged, damage reported as reason',
         'c: Previous damage, damage not reported as reason',
         'd: No damage                                   ')
  initial_status=data.frame(
    Proportion_of_Attrition=c(0,a_c/a_1,0,1-a_c/a_1),
    Attrition=c(0,a_c,0,a_1-a_c),
    group=cats,
    Age.f=rep(0,4)
  )
  toPlot1=rbind(initial_status, data.frame(
    Proportion_of_Attrition=jagsResults[1:32,3],
    Attrition=jagsResults[1:32,3]*rep(dataset$a,times=4),
    group=rep(cats,each=8),
    Age.f=rep(dataset$Age.f,times=4)))
  damagecolors = c("#FF0000FF","#0000FFFF","#FF000055","#0000FF55")
  library(stringr)
  # stacked area chart
  plot1 = ggplot(data=toPlot1, aes(x=Age.f, y=Proportion_of_Attrition, fill=group)) +
    geom_area(position = position_stack(reverse = TRUE))+
    scale_fill_manual(values=damagecolors)+
    scale_x_continuous(name = 'Months since distribution',limits=c(-2,50),breaks = c(0,6,12,18,24,30,36,42,48)) +
    scale_y_continuous(name = 'Proportion of attrited nets', limits=c(0,1.01))+
    theme_bw(base_size = 8)+ theme(legend.position = "none")
  plot2 = ggplot(data=toPlot1, aes(x=Age.f, y=Attrition, fill=group)) +
    geom_area(position = position_stack(reverse = TRUE))+
    scale_fill_manual(values=damagecolors,labels=c(
      str_wrap('Previous damage, damage reported as reason',25),
      str_wrap('Previously undamaged, damage reported as reason',25),
      str_wrap('Previous damage, damage not reported as reason',25),
      str_wrap('No damage                                   ',25)),
      guide = guide_legend(reverse = TRUE))+
    scale_x_continuous(name = 'Months since distribution',limits=c(-2,50),breaks = c(0,6,12,18,24,30,36,42,48)) +
    scale_y_continuous(name = 'Attrition rate (per year)')+
    theme_bw(base_size = 8)+ theme(legend.position = "right")+
    labs(fill="Type of Attrition")
  plot3 = ggplot(data= dataset, aes(x=Age.f, y = prop))+
    geom_point(stat = "identity",size=2) + geom_errorbar(aes(ymin = lowerL, ymax = upperL),width = 2, position=position_dodge(.9))+
    geom_line(aes(y = p_x),color='blue',size=1)+
    geom_ribbon(aes(ymin=p_x_lower,ymax=p_x_upper),fill='blue',alpha=0.3)+
    scale_x_continuous(name = 'Months since distribution',limits=c(-2,50),breaks = c(0,6,12,18,24,30,36,42,48)) +
    scale_y_continuous(name = 'Proportion of attrited nets', limits=c(0,1))+
    theme_bw(base_size = 8)
  library(cowplot)
  plot =plot_grid(
    plotlist = list(plot3,plot1,plot2),
    ncol = 3,
    rel_widths = c(1,1,1.75),
    labels="AUTO")
  savePlot(plot,'AnalysisOfRecalls.png',vertical_panels=0.6)
return(list(estimates=jagsResults,plot1=plot1,plot2=plot2,plot3=plot3))}

categories <- c("New","Used last night","Not used last night","Not in use","Attrited")
extended_categories <- c("New","Used last night","Not used last night","Not in use","Attrited","Removed","Absent")
hole_categories <- c('New','Undamaged-Not in use','Undamaged-In use','Damaged-In use','Damaged-Not in use','Attrited')
folderForDurationFiles = ''



