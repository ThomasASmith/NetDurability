CHOOSING THE RIGHT MODEL:

prior: lognormal or uniform; and each time check different parameters 
--> With uniform priors, a1 and h1 were not sufficiently 'away from zero'. The other parametes were not different for lognormal(0,2) vs. uniform(0,10) (see pairs 'netdurability_fixq_uniformall__pairs.png').

newnets model: free q or fixed q 
--> We found that the histogram of the q posterior is not sufficiently 'away from zero' (pairs plot: 'pairsplot_freeq.png'). Hence, we fix q to 0.

repair model: I currently set the counts corresponding to repaired nets to 0.

parameterisation: plain, ratio or function
--> A reparameterisation with u1, u1/(u1+v3), u2, u2/(u2+v4) instead of u1,v3,u2,v4 lead to the same posteriors of the rates. 
--> A parameterisation with u1 = b1 + alpha1*v3 and u2 = b2 + alpha2*v4 found the same posteriors for u1,u2,v3,v4 ('see netdurabiity_fixq_reparam_func__pairs.png'), but no clear estimate for alpha1,alpha2,b1,b2 ('see netdurabiity_fixq_reparam_func__pairs2.png'). Hence, no linear relation between u1 and v3 as well as u2 and v4 can be proven.


REMARKS: 
I suggest to go for the model 'netdurability'. 
But it needs to be separately shown that this model fits well the data, in particular that the Markov assumption is valid for this data.

No hierarchy. So, NODATA fit not important.

I sometimes got warnings on maximum treedepth, BFMI and ESS. I think these all came from hitting the maximal tree depth in many cases (~350). I use 4 chains with each 2.000 warmup and 2.000 post-warmup draws, otherwise I haven't taken any action. Note that the model with free q was run with 1.000 post-warmup iterations only.

My opinion on parameter correlation: The parameter pairs u1 and v3 as well as u2 and v4 are clearly positively correlated. But even though the dependency looks linear, no linear relation can be proven. And, the credible region is a narrow ellipsis but with limited spread. Therefore, I think the inference is OK since we can clearly localise all parameters and only within small regions same parameters are correlated. 

My opinion all parameters close to 0: I think the inference is OK, but a1 and h1, and even q, are simply not significant. This is a valid conclusion, suggesting that new nets are first stored and then gradually used, and that stored nets don't get holed or thrown away.

The parameter p0 represents the probability that a net that is checked for use is also checked for holes. Hence, this is an auxiliary parameter needed to analyse the data. It can be shown that the maximul likelihood estimate of the other, important, parameters is independent of p0.


TODO ONCE METHOD FOUND:
save sample, at least in list and then as csv
save print
save pairs plot
I don't plot the priors and posteriors since I don't know what format Tom wants and since I am usually doing this in matlab. 
if needed: plot chains
