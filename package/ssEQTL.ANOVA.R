# created on June 21, 2020
#  (1) write a wrapper function
#  (2) rename input parameters
#  (3) remove powerEQTL.ANOVA2 related functions
#
# created on Dec. 11, 2016
# (1) allow mu2-mu1 not equal to mu3-mu2,
#   where mu1 is the mean expression level for mutation homozygote,
# mu2 is the mean expression level for heterozygote,
#  and mu3 is the mean expression level for wildtype homozygote
#
# created on Dec. 8, 2016
#  (1) sample size calculation for eQTL analysis based on ANOVA or simple linear regression
#

# MAF - minor allele frequency
# FWER - family-wise type I error rate
# nTests - number of tests
# n - total number of subjects
# sigma - standard deviation of gene expression levels
#            (assume each group of subjects has the same mystddev)
# deltaVec - mean difference of gene expression levels between groups
#            deltaVec[1]=mu1-mu2 and deltaVec[2]=mu3-mu2
#  group 1 is mutation homozygotes
#  group 2 is heterozygotes
#  group 3 is wildtype homozygotes

# squared difference between estimated power and desired power
diffPower4ss.ANOVA=function(n,
                       MAF=0.1,
                       deltaVec=c(-0.13,0.13),
                       power=0.8,
                       sigma=0.13,
                       FWER=0.05,
                       nTests=200000)
{
  est.power=powerEQTL.ANOVA.default(MAF=MAF,
                              n=n,
                              sigma=sigma,
                              deltaVec=deltaVec, 
                              FWER=FWER,
                              nTests=nTests)
  diff=est.power - power
  return(diff)

}

# MAF - minor allele frequency
# typeI - type I error rate
# nTests - number of tests
# mypower - desired power
# sigma - standard deviation of gene expression levels
#            (assume each group of subjects has the same mystddev)
# deltaVec - mean difference of gene expression levels between groups
#            deltaVec[1]=mu2-mu1 and deltaVec[2]=mu3-mu2
ssEQTL.ANOVA=function(MAF,
                      deltaVec=c(-0.13,0.13),
                      power=0.8,
                      sigma=0.13,
                      FWER=0.05,
                      nTests=200000,
                      n.lower = 2.01,
                      n.upper = 1e+30)
{

  res.root=uniroot(f=diffPower4ss.ANOVA,
                   interval = c(n.lower, n.upper),
                   MAF=MAF,
                   deltaVec=deltaVec,
                   power=power,
                   sigma=sigma,
                   FWER=FWER,
                   nTests=nTests
                  )

  return(res.root$root)
}



