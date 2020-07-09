# created on June 22, 2020

# define function
# We assume the following linear mixed effects model to characterize
#   the association between genotype and gene expression:
#  y_{ij} = beta_{0i} + beta_1 * x_i + epsilon_{ij},
#    beta_{0i} ~ N(beta_0, sigma^2_{\beta})
#    epsilon_{ij} ~ N(0, sigma^2_{epsilon})
#
# slope - slope under alternative hypothesis
# n - number of subjects
# m - number of cells per subject
# sigma.y - standard deviation of the gene expression
# MAF - minor allele frequency (between 0 and 0.5)
# rho - intra-class correlation (i.e., correlation between y_{ij} and y_{ik})
#        rho = sigma^2_{beta} / (sigma^2_{beta}+sigma^2_{epsilon})
# FWER - family-wise type I error rate
# nTests = number of genes * number of SNPs


# difference between estimated power and desired power
diffPower4ss.scRNAseq=function( n, 
                                slope, 
                                 m,
                                power,
                                 sigma.y, 
                                 MAF=0.2, 
                                 rho=0.8, 
                                 FWER=0.05,
                                 nTests=1)
{
  est.power=powerEQTL.scRNAseq.default(
    slope = slope, 
    n = n, 
    m = m, 
    sigma.y = sigma.y, 
    MAF=MAF, 
    rho=rho, 
    FWER=FWER,
    nTests=nTests)
  diff=est.power - power
  return(diff)

}


ssEQTL.scRNAseq=function(  slope, 
                           m,
                           power = 0.8,
                           sigma.y = 0.29, 
                           MAF=0.2, 
                           rho=0.8, 
                           FWER=0.05,
                           nTests=1,
                           n.lower=2.01,
                           n.upper=1e+30)
{

  res.root=uniroot(f=diffPower4ss.scRNAseq,
                   interval = c(n.lower, n.upper),
                   slope = slope, 
                   m = m,
                   power = power,
                   sigma.y = sigma.y, 
                   MAF=MAF, 
                   rho=rho, 
                   FWER=FWER,
                   nTests=nTests
                  )

  return(res.root$root)
}



