# modified on June 21, 2020
#  (1) use exact power calculation formula
#
# created on Dec. 8, 2016
#  (1) sample size calculation for eQTL analysis based on simple linear regression
#

# MAF - minor allele frequency
# FWER - family-wise type I error rate
# nTests - number of tests
# n - total number of subjects
# power - desired power
# sigma.y - standard deviation of the outcome
#            
# n.lower - lower bound for sample size
# n.upper - upper bound for sample size
# verbose - flag indicating if intermedaite results should be output

diffPower4ss.SLR=function(n,
                          MAF,
                          slope,   
                          sigma.y=0.13,
                          FWER=0.05,
                          nTests=200000,
                          power=0.8)
{
  power.est = powerEQTL.SLR.default(MAF = MAF,
                                     slope=slope,
                                     n=n,
                                     sigma.y=sigma.y,
                                     FWER=FWER,
                                     nTests=nTests)
  diff=power.est-power
  return(diff)
  
}


ssEQTL.SLR=function(MAF,
                       slope=0.13,
                       power=0.8,
                       sigma.y=0.13,
                       FWER=0.05,
                       nTests=200000,
                       n.lower = 2.01,
                       n.upper = 1e+30
                       )
{
  sigma2.x = 2*MAF*(1-MAF)
  
  delta = slope
  bound = sigma.y/sqrt(sigma2.x)
  
  if(delta >= bound || delta <= - bound)
  {
    stop("slope must be in the interval (-a, a), where a = sigma.y/sqrt(2MAF(1-MAF))!\n")
  }
  
  res.uni=uniroot(f=diffPower4ss.SLR,
                  interval = c(n.lower, n.upper),
                  MAF = MAF,
                  slope = slope,   
                  sigma.y=sigma.y,
                  FWER=FWER,
                  nTests=nTests,
                  power=power
                  )
  return(res.uni$root)
}


