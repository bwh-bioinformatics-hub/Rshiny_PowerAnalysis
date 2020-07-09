
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


powerEQTL.ANOVA=function(MAF,
                         deltaVec=c(-0.13, 0.13),
                         n=200,
                         power = NULL,
                         sigma=0.13,
                         FWER=0.05,
                         nTests=200000,
                         n.lower = 4,
                         n.upper = 1e+30)
{
    if(is.null(MAF)==TRUE && 
       is.null(n) == FALSE && is.null(power) == FALSE)
    {
        MAF = minMAFeQTL.ANOVA(
            deltaVec=deltaVec,
            n=n,
            power = power,
            sigma=sigma,
            FWER=FWER,
            nTests=nTests)
        
        names(MAF) = "MAF"
        return(MAF)
    } else if (is.null(MAF)==FALSE &&
               is.null(n) == TRUE && is.null(power) == FALSE) {
        n = ssEQTL.ANOVA(MAF = MAF,
                         deltaVec=deltaVec,
                         power=power,
                         sigma=sigma,
                         FWER=FWER,
                         nTests=nTests,
                         n.lower = n.lower,
                         n.upper = n.upper)
        
        names(n)="n"
        return(n)
    } else if (is.null(MAF)==FALSE && 
               is.null(n) == FALSE && is.null(power) == TRUE) {
        power = powerEQTL.ANOVA.default(MAF = MAF,
                                        deltaVec=deltaVec,
                                        n=n,
                                        sigma=sigma,
                                        FWER=FWER,
                                        nTests=nTests) 
        names(power) = "power"
        return(power)
    } else {
        stop("One and only one of the 3 parameters (MAF, n, power) can be NULL!\n")
    }
    
}


powerEQTL.ANOVA.default=function(MAF,
                                 deltaVec=c(-0.13, 0.13),
                                 n=200,
                                 sigma=0.13,
                                 FWER=0.05,
                                 nTests=200000)
{
    if(length(deltaVec)!=2)
    {
        stop("'deltaVec' has 2 and only 2 elements!\n1st element = mu2 - mu1; 2nd element = mu3 - mu2!\n")
    }
    
    gm1 = -deltaVec[1] # mu2 - mu1
    gm2 = 0
    gm3 = deltaVec[2] # mu3 - mu2
    
    w1=MAF^2 # mutation homozygotes
    w2=2*MAF*(1-MAF) # heterozygotes
    w3=(1-MAF)^2 # wildtype homozygotes
    
    
    alpha=FWER/nTests
    
    k=3
    mydf1=k-1
    mydf2=n-k
    q=qf(p=1-alpha, df1=mydf1, df2=mydf2)
    
    wVec=c(w1, w2, w3)
    muVec=c(gm1, gm2, gm3)
    
    mu=sum(wVec*muVec, na.rm=TRUE)
    
    myncp = n*sum(wVec*(muVec-mu)^2, na.rm=TRUE)
    myncp=myncp/(sigma^2)
    
    power=1-pf(q=q, df1=mydf1, df2=mydf2,
               ncp=myncp)
    
    return(power)
    
}


# MAF - minor allele frequency
# FWER - family-wise type I error rate
# nTests - number of tests
# n - total number of subjects
# sigma.y - standard deviation of the outcome


powerEQTL.SLR=function(MAF,
                       slope=0.13,
                       n=200,
                       power = NULL,
                       sigma.y=0.13,
                       FWER=0.05,
                       nTests=200000,
                       n.lower = 2.01,
                       n.upper = 1e+30)
{
    if(is.null(MAF)==TRUE && is.null(slope) == FALSE &&
       is.null(n) == FALSE && is.null(power) == FALSE)
    {
        MAF = minMAFeQTL.SLR(slope=slope,
                             n=n,
                             power=power,
                             sigma.y=sigma.y,
                             FWER=FWER,
                             nTests=nTests)
        names(MAF) = "MAF"
        return(MAF)
    } else if (is.null(MAF)==FALSE && is.null(slope) == TRUE &&
               is.null(n) == FALSE && is.null(power) == FALSE) {
        slope = minSlopeEQTL.SLR(MAF = MAF,
                                 n= n,
                                 power=power,
                                 sigma.y=sigma.y,
                                 FWER=FWER,
                                 nTests=nTests)
        names(slope) = "slope"
        return(slope)
    } else if (is.null(MAF)==FALSE && is.null(slope) == FALSE &&
               is.null(n) == TRUE && is.null(power) == FALSE) {
        n = ssEQTL.SLR(MAF = MAF,
                       slope=slope,
                       power=power,
                       sigma.y=sigma.y,
                       FWER=FWER,
                       nTests=nTests,
                       n.lower = n.lower,
                       n.upper = n.upper)
        names(n)="n"
        return(n)
    } else if (is.null(MAF)==FALSE && is.null(slope) == FALSE &&
               is.null(n) == FALSE && is.null(power) == TRUE) {
        power = powerEQTL.SLR.default(MAF = MAF,
                                      slope=slope,
                                      n=n,
                                      sigma.y=sigma.y,
                                      FWER=FWER,
                                      nTests=nTests) 
        names(power) = "power"
        return(power)
    } else {
        stop("One and only one of the 4 parameters (MAF, slope, n, power) can be NULL!\n")
    }
    
}


powerEQTL.SLR.default=function(MAF,
                               slope=0.13,
                               n=200,
                               sigma.y=0.13,
                               FWER=0.05,
                               nTests=200000)
{
    sigma2.x = 2*MAF*(1-MAF)
    sigma2.y = sigma.y^2
    alpha = FWER/nTests
    
    delta = slope
    
    bound = sigma.y/sqrt(sigma2.x)
    
    if(delta >= bound || delta <= - bound)
    {
        stop("slope must be in the interval (-a, a), where a = sigma.y/sqrt(2MAF(1-MAF))!\n")
    }
    
    numer.ncp = delta *sqrt((n-1)*sigma2.x)
    denom.ncp= sqrt(sigma2.y - delta^2*sigma2.x)
    
    lambda = numer.ncp / denom.ncp
    
    mydf = n - 2
    cutoff = qt(1-alpha/2, df=mydf, ncp=0)
    
    power = 1 - pt(cutoff, df=mydf, ncp=lambda)
    power = power + pt(-cutoff, df=mydf, ncp=lambda)
    
    return(power)
    
}

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
power.eQTL.scRNAseq=function(delta, n, m, sigma.y, theta=0.2, rho=0.8, alpha=0.05)
{
    za2=qnorm(1-alpha/2)
    sigma.x=sqrt(2*theta*(1-theta))
    part0=sigma.x*delta*sqrt(m*(n-1))/(sigma.y*sqrt(1+(m-1)*rho))
    part1 = za2-part0
    
    part2 = -za2-part0
    
    power = 1- pnorm(part1) + pnorm(part2)
    
    return(power)
    
}

powerEQTL.scRNAseq=function(
    slope, 
    n, 
    m, 
    power = NULL,
    sigma.y, 
    MAF=0.2, 
    rho=0.8, 
    FWER=0.05,
    nTests=1,
    n.lower=2.01,
    n.upper=1e+30)
{
    if(is.null(MAF)==TRUE && is.null(slope) == FALSE &&
       is.null(n) == FALSE && is.null(power) == FALSE)
    {
        MAF = minMAFEQTL.scRNAseq(
            slope = slope,
            n = n, 
            m = m,
            power = power,
            sigma.y = sigma.y, 
            rho=rho, 
            FWER=FWER,
            nTests=nTests)
        
        names(MAF) = "MAF"
        return(MAF)
    } else if (is.null(MAF)==FALSE && is.null(slope) == TRUE &&
               is.null(n) == FALSE && is.null(power) == FALSE) {
        slope = minSlopeEQTL.scRNAseq(
            n = n, 
            m = m,
            power = power,
            sigma.y = sigma.y, 
            MAF=MAF, 
            rho=rho, 
            FWER=FWER,
            nTests=nTests)
        
        names(slope) = "slope"
        return(slope)
    } else if (is.null(MAF)==FALSE && is.null(slope) == FALSE &&
               is.null(n) == TRUE && is.null(power) == FALSE) {
        n = ssEQTL.scRNAseq(  
            slope = slope, 
            m = m,
            power = power,
            sigma.y = sigma.y, 
            MAF=MAF, 
            rho=rho, 
            FWER=FWER,
            nTests=nTests,
            n.lower=n.lower,
            n.upper=n.upper)
        
        names(n)="n"
        return(n)
    } else if (is.null(MAF)==FALSE && is.null(slope) == FALSE &&
               is.null(n) == FALSE && is.null(power) == TRUE) {
        power = powerEQTL.scRNAseq.default(
            slope = slope, 
            n = n, 
            m = m, 
            sigma.y = sigma.y, 
            MAF=MAF, 
            rho=rho, 
            FWER=FWER,
            nTests=nTests)
        
        names(power) = "power"
        return(power)
    } else {
        stop("One and only one of the 4 parameters (MAF, slope, n, power) can be NULL!\n")
    }
    
}


powerEQTL.scRNAseq.default=function(
    slope, 
    n, 
    m, 
    sigma.y, 
    MAF=0.2, 
    rho=0.8, 
    FWER=0.05,
    nTests=1)
{
    alpha2=FWER/nTests
    za2=qnorm(1-alpha2/2)
    
    sigma.x=sqrt(2*MAF*(1-MAF))
    part0=sigma.x*slope*sqrt(m*(n-1))/(sigma.y*sqrt(1+(m-1)*rho))
    part1 = za2-part0
    
    part2 = -za2-part0
    
    power = 1- pnorm(part1) + pnorm(part2)
    
    return(power)
    
}
