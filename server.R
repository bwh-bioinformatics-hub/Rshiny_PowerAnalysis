library(shiny)
library(devtools)
library('powerEQTL') 
library(dplyr)
library(ggplot2)


powerEQTL.ANOVA=function(MAF, 
                   alpha=NULL,
                   n=200,
                   sigma=0.13,
                   delta=0.13)
{
    gm1 = delta # mu2 - mu1
    gm2 = 0
    gm3 = delta # mu3 - mu2
    
    w1=MAF^2 # mutation homozygotes
    w2=2*MAF*(1-MAF) # heterozygotes
    w3=(1-MAF)^2 # wildtype homozygotes
    
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
    
    lambda.a = slope
    
    numer.xi = lambda.a *sqrt((n-1)*sigma2.x)
    denom.xi= sqrt(abs(sigma2.y - lambda.a^2*sigma2.x))
    
    xi = numer.xi / denom.xi
    
    mydf = n - 2
    cutoff = qt(1-alpha/2, df=mydf, ncp=0)
    
    power = 1 - pt(cutoff, df=mydf, ncp=xi)
    power = power + pt(-cutoff, df=mydf, ncp=xi)
    
    return(power)
    
}


diffPower4ss.ANOVA=function(n,
                            MAF=0.1,
                            delta=0.2,
                            power=0.8,
                            sigma=0.2,
                            FWER=0.05,
                            nTests=200000)
{
    est.power=powerEQTL.ANOVA(MAF=MAF,
                              n=n,
                              sigma=sigma,
                              delta=delta, 
                              alpha=FWER/nTests)
    diff=est.power - power
    return(diff)
    
}

diffPower4ss.SLR=function(n,
                            MAF=0.1,
                            slope=0.13,
                            power=0.8,
                            sigma=0.13,
                            FWER=0.05,
                            nTests=200000)
{
    est.power=powerEQTL.SLR.default(MAF = MAF,
                                    slope=slope,
                                    n=n,
                                    sigma.y=sigma,
                                    FWER=FWER,
                                    nTests=nTests) 
    diff=est.power - power
    return(diff)
    
}

diffPower4ss.scRNAseq=function(n,
                               m,
                          MAF=0.1,
                          slope=0.13,
                          power=0.8,
                          sigma=0.13,
                          FWER=0.05,
                          nTests=200000,
                          rho=0.8)
{
    est.power=power.eQTL.scRNAseq(theta = MAF,
                                    delta=slope,
                                    n=n,
                                    m=m,
                                    sigma.y=sigma,
                                    alpha = FWER/nTests,
                                    rho=rho) 
    diff=est.power - power
    return(diff)
    
}


function(input, output) {

    ## values for theta (MAF) (single-cell)
    mafVec=seq(from=0.005,to=0.5,by=0.001)
    mafLen=length(mafVec)

    # values for theta (MAF) (tissue)
    mafVec1 = seq(from=0.005, to=0.5, by=0.001)
    mafLen1 = length(mafVec1)
    
    ## values for sample sizes (tissue eqtl)
    subVec.data <- reactive({unlist(strsplit(as.character(input$subjects), ","))})
    sizeVec.data <- reactive({as.numeric(unlist(strsplit(as.character(input$myntot), ",")))})
    sizeLen.data <- reactive(length(sizeVec.data()))
    
    # values for sample sizes (single-cell eqtl)
    subVec1.data <- reactive({unlist(strsplit(as.character(input$subjects1), ","))})
    sizeVec1.data <- reactive({as.numeric(unlist(strsplit(as.character(input$myntot1), ",")))})
    sizeLen1.data <- reactive(length(sizeVec1.data()))
    
    ## validate intra-class correlation
    text.data <- reactive({
         shiny::validate(need((input$rho <= 1 && input$rho >= 0),
                       "Intra-class correlation must in range [0, 1]"))
         input$rho
    })
    
    text1.data <- reactive({
        shiny::validate(need((input$rho2 <= 1 && input$rho2 >= 0),
                             "Intra-class correlation must in range [0, 1]"))
        input$rho2
    })
    
    
    ## single-cell eQTL
    sc.data <- reactive({
        power_sc_eQTL <- array(numeric(sizeLen1.data()*mafLen), dim=c(sizeLen1.data(),mafLen), dimnames = list(as.character(sizeVec1.data()), as.character(mafVec)))
        for(i in 1:sizeLen1.data()) {
            for(j in 1:mafLen) {
                power_sc_eQTL[i,j]=power.eQTL.scRNAseq(delta=input$delta, n=sizeVec1.data()[i], m=input$m,sigma.y=input$sigma, theta=mafVec[j], rho=text.data(), alpha=input$alpha/input$nTest)
            }
        }
        
        xrange <- range(mafVec*100)
        yrange <- c(0:1)
        colors <- hcl.colors(sizeLen.data(),"Set 2")
        title <- "Single-cell eQTL"
        
        
        # treating each cell independantly 
        plot(xrange, yrange, xlim=range(xrange), log='x', type="n",
             xlab="MAF (%)",
             ylab="Power",
             main= title)
        
        mtext(paste("linear mixed effects model (alpha=",input$alpha/input$nTest,", eQTL effect size=",input$delta,", n_cell="
                    , input$m ,", sigma=", input$sigma ,", rho=", input$rho ,")", sep = ""),3)
        
        abline(v=0, h=seq(0,1,.1), lty=2, col="grey89")
        abline(h=0, v=c(1:10), lty=2,col="grey89")
        abline(h=0.8, col=2)
        
        # add power curves
        for (i in 1:sizeLen1.data()){
            lines(mafVec*100, power_sc_eQTL[i,], type="l", lwd=4, col=colors[i])
            pos=which(power_sc_eQTL[i,]>=0.8)[1]
            abline(v=mafVec[pos]*100, col=colors[i])
            text(mafVec[pos]*100,0.8,labels = as.character(mafVec[pos]*100), srt = 90,cex=0.6, adj=0.5)
        }
        legend("bottomright", title="Sample size (n)", paste(subVec.data(), " (n = ",sizeVec.data(),")", sep=""),
               title.col='black',text.col=colors,cex =.7, bty='n')
    })
     
    ## plot single-cell eQTL
    output$cell <- renderPlot({
        sc.data()
    })
    
    output$cell_info <- renderText({
        # nearPoints(sc.data, input$plot_click) %>%
        #     transmute(
        #         probe,
        #         gene = HUGO.gene.symbol,
        #         `eQTL effect size = ` = signif(deltaVec, digits = 2),
        #         `Power(%) = ` = signif(power_sc_eQTL, digits = 2)
        #     )
        paste0("eQTL effect size = ", input$cpoint$x, "\nPower(%) = ", input$cpoint$y)
    })
    
    ## tissue eQTL
    ## plot tissue eQTL
    output$tissue <- renderPlot({
        power_unbalanced <- array(numeric(sizeLen.data()*mafLen1), dim=c(sizeLen.data(),mafLen1), dimnames = list(sizeVec.data(), as.character(mafVec1*100)))
        
        # ANOVA 
        if (input$radio == "One-way unbalanced anova")
        {
            for (i in 1:sizeLen.data()){
                for (j in 1:mafLen1){
                    # treating each cell independantly 
                    power_unbalanced[i,j] <- powerEQTL.ANOVA(MAF=mafVec1[j], alpha=input$alpha1/input$nTest1, n=sizeVec.data()[i], delta = input$delta1, sigma = input$sigma1)
                }
            }
            
            ## plot
            
            # set up graph
            xrange <- range(mafVec1*100)
            yrange <- c(0:1)
            colors <- hcl.colors(sizeLen.data(),"Set 2")
            title <- "Tissue eQTL"
            
            ## 	 # treating each cell independantly 
            plot(xrange, yrange, xlim=range(xrange), log='x', type="n",
                 xlab="MAF (%)",
                 ylab="Power",
                 main= title)
            mtext(paste("eQTL p-value = ", input$alpha1/input$nTest1, " (one-way unbalanced ANOVA)", sep = ""),3)
            
        } else{
            # SLR
            for (i in 1:sizeLen.data()){
                for (j in 1:mafLen1){
                    # treating each cell independantly 
                    power_unbalanced[i,j] <- powerEQTL.SLR.default(MAF=mafVec1[j], FWER=input$alpha1, nTests = input$nTest1, n=sizeVec.data()[i], slope = input$delta1, sigma.y = input$sigma1)
                }
            }
            
            ## plot
            
            # set up graph
            xrange <- range(mafVec1*100)
            yrange <- c(0:1)
            colors <- hcl.colors(sizeLen.data(),"Set 2")
            title <- "Tissue eQTL"
            
            ## 	 # treating each cell independantly 
            plot(xrange, yrange, xlim=range(xrange), log='x', type="n",
                 xlab="MAF (%)",
                 ylab="Power",
                 main= title)
            mtext(paste("eQTL p-value = ", input$alpha1/input$nTest1, " (Simple Linear Regression)", sep = ""),3)
        }
        
        abline(v=0, h=seq(0,1,.1), lty=2, col="grey89")
        abline(h=0, v=c(1:10), lty=2,col="grey89")
        abline(h=0.8, col=2)
        # add power curves
        for (i in 1:sizeLen.data()){
            lines(mafVec1*100, power_unbalanced[i,], type="l", lwd=4, col=colors[i])
            pos1=which(power_unbalanced[i,]>=0.8)[1]
            abline(v=mafVec1[pos1]*100, col=colors[i])
            text(mafVec1[pos1]*100,0.8,labels = as.character(mafVec1[pos1]*100), srt = 90,cex=0.6, adj=0.5)
        }
        legend("bottomright", title="Sample size (n)", paste(subVec.data(), " (n = ",sizeVec.data(),")", sep=""),
               title.col='black',text.col=colors,cex =.7, bty='n')
    })
    
    output$tissue_info <- renderText({
        paste0("MAF(%) = ", input$plot_click1$x, "\nPower(%) = ", input$plot_click1$y, sep = "")
    })
    
    ## clicking on the export button will generate a pdf file for sc eQTL
    output$export = downloadHandler(
        filename = function() {"single_cell_eqtl.pdf"},
        content = function(file) {
            pdf(file, width = 10, height = 8)
            
            power_sc_eQTL <- array(numeric(sizeLen1.data()*mafLen), dim=c(sizeLen1.data(),mafLen), dimnames = list(as.character(sizeVec1.data()), as.character(mafVec)))
            for(i in 1:sizeLen1.data()) {
                for(j in 1:mafLen) {
                    power_sc_eQTL[i,j]=power.eQTL.scRNAseq(delta=input$delta, n=sizeVec1.data()[i], m=input$m,sigma.y=input$sigma, theta=mafVec[j], rho=text.data(), alpha=input$alpha)
                }
            }
            
            xrange <- range(mafVec*100)
            yrange <- c(0:1)
            colors <- hcl.colors(sizeLen.data(),"Set 2")
            title <- "Single-cell eQTL"
            
            
            # treating each cell independantly 
            plot(xrange, yrange, xlim=range(xrange), log='x', type="n",
                 xlab="MAF (%)",
                 ylab="Power",
                 main= title)
            
            mtext(paste("linear mixed effects model (alpha=",input$alpha,", eQTL effect size=",input$delta,", n_cell="
                        , input$m ,", sigma=", input$sigma ,", rho=", input$rho ,")", sep = ""),3)
            
            abline(v=0, h=seq(0,1,.1), lty=2, col="grey89")
            abline(h=0, v=c(1:10), lty=2,col="grey89")
            abline(h=0.8, col=2)
            
            # add power curves
            for (i in 1:sizeLen1.data()){
                lines(mafVec*100, power_sc_eQTL[i,], type="l", lwd=4, col=colors[i])
                pos=which(power_sc_eQTL[i,]>=0.8)[1]
                abline(v=mafVec[pos]*100, col=colors[i])
                text(mafVec[pos]*100,0.8,labels = as.character(mafVec[pos]*100), srt = 90,cex=0.6, adj=0.5)
            }
            legend("bottomright", title="Sample size (n)", paste(subVec.data(), " (n = ",sizeVec.data(),")", sep=""),
                   title.col='black',text.col=colors,cex =.7, bty='n')
            
            
            dev.off()
        }
    )
    
    ## clicking on the export button will generate a pdf file for tissue eQTL
    output$export2 = downloadHandler(
        filename = function() {"tissue_eqtl.pdf"},
        content = function(file) {
            pdf(file, width = 10, height = 8)
            
            power_unbalanced <- array(numeric(sizeLen.data()*mafLen1), dim=c(sizeLen.data(),mafLen1), dimnames = list(sizeVec.data(), as.character(mafVec1*100)))
            
            # ANOVA 
            if (input$radio == "One-way unbalanced anova")
            {
                for (i in 1:sizeLen.data()){
                    for (j in 1:mafLen1){
                        # treating each cell independantly 
                        power_unbalanced[i,j] <- powerEQTL.ANOVA(MAF=mafVec1[j], alpha=input$alpha1/input$nTest1, n=sizeVec.data()[i], delta = input$delta1, sigma = input$sigma1)
                    }
                }
                
                ## plot
                
                # set up graph
                xrange <- range(mafVec1*100)
                yrange <- c(0:1)
                colors <- hcl.colors(sizeLen.data(),"Set 2")
                title <- "Tissue eQTL"
                
                ## 	 # treating each cell independantly 
                plot(xrange, yrange, xlim=range(xrange), log='x', type="n",
                     xlab="MAF (%)",
                     ylab="Power",
                     main= title)
                mtext(paste("eQTL p-value = ", input$alpha1/input$nTest1, " (one-way unbalanced ANOVA)", sep = ""),3)
                
            } else{
                # SLR
                for (i in 1:sizeLen.data()){
                    for (j in 1:mafLen1){
                        # treating each cell independantly 
                        power_unbalanced[i,j] <- powerEQTL.SLR.default(MAF=mafVec1[j], FWER=input$alpha1, nTests = input$nTest1, n=sizeVec.data()[i], slope = input$delta1, sigma.y = input$sigma1)
                    }
                }
                
                ## plot
                
                # set up graph
                xrange <- range(mafVec1*100)
                yrange <- c(0:1)
                colors <- hcl.colors(sizeLen.data(),"Set 2")
                title <- "Tissue eQTL"
                
                ## 	 # treating each cell independantly 
                plot(xrange, yrange, xlim=range(xrange), log='x', type="n",
                     xlab="MAF (%)",
                     ylab="Power",
                     main= title)
                mtext(paste("eQTL p-value = ", input$alpha1/input$nTest1, " (Simple Linear Regression)", sep = ""),3)
            }
            
            abline(v=0, h=seq(0,1,.1), lty=2, col="grey89")
            abline(h=0, v=c(1:10), lty=2,col="grey89")
            abline(h=0.8, col=2)
            # add power curves
            for (i in 1:sizeLen.data()){
                lines(mafVec1*100, power_unbalanced[i,], type="l", lwd=4, col=colors[i])
                pos1=which(power_unbalanced[i,]>=0.8)[1]
                abline(v=mafVec1[pos1]*100, col=colors[i])
                text(mafVec1[pos1]*100,0.8,labels = as.character(mafVec1[pos1]*100), srt = 90,cex=0.6, adj=0.5)
            }
            legend("bottomright", title="Sample size (n)", paste(subVec.data(), " (n = ",sizeVec.data(),")", sep=""),
                   title.col='black',text.col=colors,cex =.7, bty='n')
            
            dev.off()
        }
    )
    
    ## Input Entries Explanations
    
    output$effectSize <- renderPlot({
        gene <- data.frame(x = c( rep("CC", 20), rep("TC", 20), rep("TT", 20)), 
                           y = c( rnorm(20, 0.05, 0.1), rnorm(20, 0.0, 0.05), rnorm(20, -0.1, 0.2)))
        gene %>% ggplot( aes(x = x, y = y, fill = x)) +
            geom_boxplot() +
            geom_abline(aes(intercept = 0.125, slope = -0.07)) +
            annotate("text", label = "slope = eQTL effect size", x = 1.5, y = -0.2) +
            xlab("") +
            ylab("Rank Normalized Gene Expression") +
            ggtitle("cis-eQTL analysis") +
            theme(plot.title = element_text(lineheight = 0.8, face = "bold"))
    })
    
    output$effectSize1 <- renderPlot({
        gene <- data.frame(x = c( rep("CC", 20), rep("TC", 20), rep("TT", 20)), 
                           y = c( rnorm(20, 0.05, 0.1), rnorm(20, 0.0, 0.05), rnorm(20, -0.1, 0.2)))
        gene %>% ggplot( aes(x = x, y = y, fill = x)) +
            geom_boxplot() +
            geom_abline(aes(intercept = 0.125, slope = -0.07)) +
            annotate("text", label = "slope = eQTL effect size", x = 1.5, y = -0.2) +
            xlab("") +
            ylab("Rank Normalized Gene Expression") +
            ggtitle("cis-eQTL analysis") +
            theme(plot.title = element_text(lineheight = 0.8, face = "bold"))
    })
    
    output$description1 <- renderUI({
        withMathJax(helpText("$$\\text{intra-class correlation = correlation between } y_{ij} \\text{ and } y_{ik}$$"),
                    helpText("$$\\text{= }\\sigma^2_{\\beta} / (\\sigma^2_{\\beta}+\\sigma^2_{\\epsilon})$$"),
                    # helpText("$$\\text{= # of folds * standard deviation}$$"),
                    helpText("$$\\text{alpha = family-wise type I error rate/nTests}$$"),
                    helpText("$$\\text{nTests = number of genes * number of SNPs}$$"))
    })
    
    output$description <- renderUI({
        tags$div(
        
        withMathJax(
            helpText("Power calculation for association between genotype and gene expression based on single cell RNAseq data."),
            helpText("We assume the following simple linear mixed effects model for each (SNP, gene) pair to characterize the association between genotype and gene expression:"),
            helpText(" $$y_{ij} = \\beta_{0i} + \\beta_1  x_i + \\epsilon_{ij}$$"),
            helpText(" $$\\beta_{0i} \\text{ ~ } N(\\beta_0, \\sigma^2_{\\beta})$$ "),
            helpText(" $$\\epsilon_{ij} \\text{ ~ } N(0, \\sigma^2_{\\epsilon})$$"),
            helpText("$$i = 1 ... n$$"),
            helpText("$$j = 1 ... m$$"),
            helpText("n is the number of subjects, m is the number of cells per subject."),
            helpText("More details can be found in the supplementary documents")
            ),
        
        tags$b("References"),
        tags$br(),
        
        paste0("Dong X, Xiaoqi L., and Qiu W. Power Calculation for Association Between Genotype and Gene Expression Based on Single Cell RNAseq Data. manuscript. (2020)")
        )
    })
    
    output$description3 <- renderUI({
        withMathJax(
            helpText("$$\\text{alpha = family-wise type I error rate/nTests}$$"),
            helpText("$$\\text{nTests = number of genes * number of SNPs}$$"))
    })
    
    output$description2 <- renderUI({
        tags$div(
            tags$h5("One-way Unbalanced ANOVA Model"),
            
            withMathJax(),
      
            helpText("Power calculation for eQTL analysis that tests if a SNP is associated to a gene probe by using unbalanced one-way ANOVA."),
            helpText("The assumption of the ANOVA approach is that the association of a SNP to a gene probe is tested by using un-balanced one-way 
                     ANOVA (e.g. Lonsdale et al. 2013). According to SAS online document 
                     https://support.sas.com/documentation/cdl/en/statug/63033/HTML/default/viewer.htm#statug_power_a0000000982.htm, 
                     the power calculation formula is:"),
            helpText(" $$power = Pr(F \\ge F_{1-\\alpha}(k-1, N-k) | F\\text{ ~ }F_{k-1, N-k, \\lambda})$$"),
            helpText(" $$\\lambda = \\frac{N}{\\sigma^2 } \\sum_{i = 1}^k w_i (\\mu_i - \\mu)^2$$"),
            helpText("$$\\mu = \\sum_{i = 1}^k w_i \\mu_i$$"),
            helpText("where k = 3 is the number of groups of subjects, N is the total number of subjects."),
            helpText("More details can be found in the supplementary documents."),
        
            tags$h5("Simple Linear Regression"),
            
            helpText("Power calculation for eQTL analysis that tests if a SNP is associated to a gene probe by using simple linear regression."),
            helpText("To test if a SNP is associated with a gene probe, we use the simple linear regression based on Dupont and Plummer (1998):"),
            helpText(" $$y_{i} = \\gamma + \\lambda  x_i + \\epsilon_{i}$$"),
            helpText(" $$\\epsilon_{i} \\text{ ~ } N(0, \\rho^2)$$"),
            helpText("$$i = 1 ... n$$"),
            helpText("More details can be found in the supplementary documents."),
            
            tags$b("References"),
            tags$br(),
            
            paste0("Dupont, W.D. and Plummer, W.D.. Power and Sample Size Calculations for Studies Involving Linear Regression. Controlled Clinical Trials. 1998;19:589-601."),
            tags$br(),
            paste0("Lonsdale J and Thomas J, et al. The Genotype-Tissue Expression (GTEx) project. Nature Genetics, 45:580-585, 2013.")
            )
    })
    
    
    ## Approximate sample sizes with given inputs
    # ANOVA model
    # output$radio = renderUI{
    #     radioButtons(inputId = "models", label = "Tissue eQTL - Select model", choices = c("One-way unbalanced anova", "Simple linear regression"))
    # }
    
    ## Check inputs
    maf.data <- reactive({
        shiny::validate(need((!is.na(input$maf2)),
                             "Need parameters MAF not empty"))
        input$maf2
    })
    power.data <- reactive({
        shiny::validate(need((!is.na(input$power) && input$power <= 1 && input$power >= 0),
                             "Need parameters power not emptyand in range [0, 1]"))
        input$power
    })
    eSize.data <- reactive({
        shiny::validate(need((!is.na(input$eSize)),
                             "Need parameters effect Size not empty"))
        input$eSize
    })
    m.data <- reactive({
        shiny::validate(need((!is.na(input$m3)),
                             "Need parameters m not empty"))
        input$m3
    })
    rho.data <- reactive({
        shiny::validate(need((!is.na(input$rho3) && input$rho3 <= 1 && input$rho3 >= 0),
                             "Need parameters intra-class correlation not empty and in range [0, 1]"))
        input$rho3
    })
    
    ## First start from an initial guess: n = 100
    ## Then increment n by 100 in each iteration to approximate n when power level approaches input with error < 0.001
    ## Each time power level > input, cut increment value by 1/2
    ## Power level shall be < 1.0
    output$approx <- renderTable({
        n=uniroot(f=diffPower4ss.ANOVA,
                         interval = c(10, 1e30),
                         MAF=maf.data(),
                         delta=eSize.data() * 0.2,
                         power=power.data(),
                         sigma=0.2,
                         FWER=0.05,
                         nTests=200000
        )
        
        n2=uniroot(f=diffPower4ss.SLR,
                  interval = c(10, 1e30),
                  MAF=maf.data(),
                  slope=eSize.data() * 0.2,
                  power=power.data(),
                  sigma=0.2,
                  FWER=0.05,
                  nTests=200000
        )
        
        data.frame(Model_used = c("One-way unbalanced Anova", "Simple Linear Regression"),
                   Number_of_subjects_needed = c(n$root, n2$root),
                   Minor_Allele_Frequency = rep(input$maf2,2),
                   Power_level = rep(input$power,2),
                   eQTL_effect_size = rep(input$eSize,2)
                   )
    })
    
    output$approx1 <- renderTable({
        n=uniroot(f=diffPower4ss.scRNAseq,
                  interval = c(10, 1e30),
                  MAF=maf.data(),
                  m=m.data(),
                  rho=rho.data(),
                  slope=eSize.data() * 0.2,
                  power=power.data(),
                  sigma=0.2,
                  FWER=0.05,
                  nTests=200000
        )
        
        data.frame(Model_used = c("Simple linear mixed effects model"),
                   No_subjects_needed = n$root,
                   Minor_Allele_Frequency = input$maf2,
                   Power_level = input$power,
                   eQTL_effect_size = input$eSize,
                   m = input$m3,
                   Intra_class_correlation = input$rho3
        )
    })
    
    output$title <- renderUI(tags$h5("Tissue eQTL"))
    
    output$title1 <- renderUI(tags$div(tags$br(), tags$h5("Single-cell eQTL")))
    
    output$Explanation <- renderUI({
        tags$div(
            withMathJax(),
            tags$br(),
            helpText("alpha = Family Type-I Error Rate / Number of Tests = 0.05 / 200000 = 5e^-8"),
            # tags$br(),
            helpText("eQTL effect size = mean difference of gene expression levels between groups(delta) / standard deviation of gene expression levels(sigma)"),
            tags$b("References"),
            tags$br(),
            
            paste("Dong X, Xiaoqi L., and Qiu W. Power Calculation for Association Between Genotype and Gene Expression Based on Single Cell RNAseq Data. manuscript. (2020)"),
            tags$br(),
            paste0("Dupont, W.D. and Plummer, W.D.. Power and Sample Size Calculations for Studies Involving Linear Regression. Controlled Clinical Trials. 1998;19:589-601."),
            tags$br(),
            paste0("Lonsdale J and Thomas J, et al. The Genotype-Tissue Expression (GTEx) project. Nature Genetics, 45:580-585, 2013.")
            
            )
    })
}



# 
# if (input$models == "One-way unbalanced anova")
# {
#     numericInput(inputId = "n2", label = "Sample size (n)",
#                  value = 640, min = 0)
#     numericInput(inputId = "maf2", label = "Minor allele frequencies (between 0 and 0.5)",
#                  value = 0.02)
#     numericInput(inputId = "power", label = "Desired Power level",
#                  value = 0.8, min = 0, max = 1)
#     
#     submitButton("submit")
#     
#     tableOutput("tANOVA")
# }
# if (input$models == "Simple linear regression")
# {
#     numericInput(inputId = "n2", label = "Sample size (n)",
#                  value = 640, min = 0)
#     numericInput(inputId = "maf2", label = "Minor allele frequencies (between 0 and 0.5)",
#                  value = 0.02)
#     numericInput(inputId = "power", label = "Desired Power level",
#                  value = 0.8, min = 0, max = 1)
#     numericInput(inputId = "delta2", label = "eQTL effect size", 
#                  value = 0.29*1.5)
#     
#     submitButton("submit")
#     
#     tableOutput("tSLR")
#     
# }


# 
# # n & inc & approx for single cell eqtl
# # n2 & inc2 & approx2 for tissue eqtl
# n <- 100
# n2 <- 100
# inc <- 100
# inc2 <- 100
# approx <- power.eQTL.scRNAseq(delta=input$delta2, n=n, m=input$m2,sigma.y=input$sigma2, theta=input$maf2, rho=text1.data(), alpha=input$alpha2)
# approx2 <- powerEQTL(MAF=input$maf2, alpha=input$alpha2, myntotal=n2, delta = input$delta2)
# 
# while (isTRUE(approx - input$power < 0 || approx - input$power > 0.005)) 
# {
#     if (approx > input$power)
#     {
#         n = n - inc
#         inc = inc/2
#     }
#     n = n + inc
#     approx = power.eQTL.scRNAseq(delta=input$delta2, n=n, m=input$m2,sigma.y=input$sigma2, theta=input$maf2, rho=text1.data(), alpha=input$alpha2)
# }
# 
# while (isTRUE(approx2 - input$power < 0 || approx2 - input$power > 0.005)) 
# {
#     if (approx2 > input$power)
#     {
#         n2 = n2 - inc2
#         inc2 = inc2/2
#     }
#     n2 = n2 + inc2
#     approx2 = powerEQTL(MAF=input$maf2, alpha=input$alpha2, myntotal=n2, delta = input$delta2)
# }
# 
# withMathJax(
#     helpText("The minium sample size to approximate power level ", input$delta2, "is ", ceiling(n))
#     ## the outputs for two models are the same
#     # helpText("The minium sample size to approximate power level using one-way unbalanced ANOVA model", input$delta2, " is ", ceiling(n2))
# )