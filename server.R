library(shiny)
library(devtools)
library('powerEQTL') 
library(dplyr)
library(ggplot2)


powerEQTL=function(MAF, 
                   alpha=NULL,
                   typeI=NULL, 
                   nSNPs=NULL, 
                   myntotal=200,
                   mystddev=0.2908849,
                   delta=0.2908849)
{
    # check parameters
    if(is.null(alpha)){
        if(is.null(typeI) || is.null(nSNPs)) stop("Either alpha or (typeI and nSNPs) should be specified.") else alpha=typeI/nSNPs;
    }
    #cat("alpha=", alpha, "\n")
    
    gm1=3*delta
    gm2=2*delta
    gm3=delta
    
    #cat("mu2=", gm1, ", mu1=", gm2, ", mu0=", gm3, "\n")
    
    n2=MAF^2
    n1=2*MAF*(1-MAF)
    n0=(1-MAF)^2
    
    #cat("n2=", n2, ", n1=", n1, ", n0=", n0, "\n")
    
    w2=1
    w1=n1/n2
    w0=n0/n2
    
    k=3
    mydf1=k-1
    mydf2=myntotal-k
    q=qf(p=1-alpha, df1=mydf1, df2=mydf2)
    
    wVec=c(w2, w1, w0)
    wVec=wVec/sum(wVec, na.rm=TRUE)
    muVec=c(gm1, gm2, gm3)
    
    mu=sum(wVec*muVec, na.rm=TRUE)
    
    myncp = myntotal*sum(wVec*(muVec-mu)^2, na.rm=TRUE)
    myncp=myncp/(mystddev^2)
    
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

function(input, output) {
    
    # value for delta
    deltaVec=seq(from=0.1, to=1, by=0.01)
    len=length(deltaVec)
    
    # values for theta (MAF)
    mafVec.data <- reactive({as.numeric(unlist(strsplit(as.character(input$MAF), ",")))})
    mafLen.data <- reactive(length(mafVec.data()))

    # values for theta (MAF) (tissue)
    mafVec1 = seq(from=0.005, to=0.5, by=0.001)
    mafLen1 = length(mafVec1)
    
    # values for sample sizes
    subVec.data <- reactive({unlist(strsplit(as.character(input$subjects), ","))})
    sizeVec.data <- reactive({as.numeric(unlist(strsplit(as.character(input$myntot), ",")))})
    sizeLen.data <- reactive(length(sizeVec.data()))
    
    # validate intra-class correlation
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
    
    # #observe Events
    # vals <- reactiveValues()
    # observeEvent(input$plot_click)
    
    # single-cell eQTL
    sc.data <- reactive({
        power_sc_eQTL <- array(numeric(len*mafLen.data()), dim=c(len,mafLen.data()), dimnames = list(as.character(deltaVec), as.character(mafVec.data())))
        for(i in 1:len) {
            for(j in 1:mafLen.data()) {
                power_sc_eQTL[i,j]=power.eQTL.scRNAseq(delta=deltaVec[i], n=input$n, m=input$m,sigma.y=input$sigma, theta=mafVec.data()[j], rho=text.data(), alpha=input$alpha)
            }
        }
        
        xrange <- c(0:1)
        yrange <- c(0:1)
        colors <- hcl.colors(mafLen.data(),"Set 2")
        title <- "Single-cell eQTL"
        
        plot(xrange, yrange, xlim=range(xrange), type="n",
             xlab="eQTL effect size",
             ylab="Power",
             main=title)
        
        mtext(paste("linear mixed effects model (alpha=",input$alpha,", n_subject=",input$n,", n_cell="
                    , input$m ,", sigma=", input$sigma ,", rho=", input$rho ,")", sep = ""),3)
        abline(h=0.8, col=2)
        
        # add power curves
        for (i in 1:mafLen.data()){
            lines(deltaVec, power_sc_eQTL[,i], type="l", lwd=4, col=colors[i])
            pos=which(power_sc_eQTL[,i]>=0.8)[1]
            print(max(power_sc_eQTL[,i]))
            abline(v=deltaVec[pos], col=colors[i])
            text(deltaVec[pos],0.8,labels = as.character(deltaVec[pos]), srt = 90,cex=0.6, adj=0.5)
        }
        legend("topleft", title="MAF", as.character(mafVec.data()), lty=1, col=colors, lwd=4, title.col='black',bty='n')
        
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
    
    # tissue eQTL
    ts.data <- reactive({
        power_unbalanced <- array(numeric(sizeLen.data()*mafLen1), dim=c(sizeLen.data(),mafLen1), dimnames = list(sizeVec.data(), as.character(mafVec1*100)))

        for (i in 1:sizeLen.data()){
            for (j in 1:mafLen1){
                # treating each cell independantly 
                power_unbalanced[i,j] <- powerEQTL(MAF=mafVec1[j], alpha=input$alpha1, myntotal=sizeVec.data()[i], delta = input$delta)
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
        mtext(paste("eQTL p-value = ", input$alpha1, " (one-way unbalanced ANOVA)", sep = ""),3)
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
    
    ## plot tissue eQTL
    output$tissue <- renderPlot({
        ts.data()
    })
    
    output$tissue_info <- renderText({
        paste0("MAF(%) = ", input$plot_click1$x, "\nPower(%) = ", input$plot_click1$y, sep = "")
    })
    
    ## clicking on the export button will generate a pdf file for sc eQTL
    output$export = downloadHandler(
        filename = function() {"single_cell_eqtl.pdf"},
        content = function(file) {
            pdf(file, width = 10, height = 8)
            
            power_sc_eQTL <- array(numeric(len*mafLen.data()), dim=c(len,mafLen.data()), dimnames = list(as.character(deltaVec), as.character(mafVec.data())))
            for(i in 1:len) {
                for(j in 1:mafLen.data()) {
                    power_sc_eQTL[i,j]=power.eQTL.scRNAseq(delta=deltaVec[i], n=input$n, m=input$m,sigma.y=input$sigma, theta=mafVec.data()[j], rho=text.data(), alpha=input$alpha)
                }
            }
            
            xrange <- c(0:1)
            yrange <- c(0:1)
            colors <- hcl.colors(mafLen.data(),"Set 2")
            title <- "Single-cell eQTL"
            
            plot(xrange, yrange, xlim=range(xrange), type="n",
                 xlab="Slope",
                 ylab="Power",
                 main=title)
            
            mtext(paste("linear mixed effects model (alpha=",input$alpha,", n_subject=",input$n,", n_cell="
                        , input$m ,", sigma=", input$sigma ,", rho=", input$rho ,")", sep = ""),3)
            abline(h=0.8, col=2)
            
            # add power curves
            for (i in 1:mafLen.data()){
                lines(deltaVec, power_sc_eQTL[,i], type="l", lwd=4, col=colors[i])
                pos=which(power_sc_eQTL[,i]>=0.8)[1]
                abline(v=deltaVec[pos], col=colors[i])
                text(deltaVec[pos],0.8,labels = as.character(deltaVec[pos]), srt = 90,cex=0.6, adj=0.5)
            }
            legend("topleft", title="MAF", as.character(mafVec.data()), lty=1, col=colors, lwd=4, title.col='black',bty='n')
            
            
            dev.off()
        }
    )
    
    ## clicking on the export button will generate a pdf file for tissue eQTL
    output$export2 = downloadHandler(
        filename = function() {"tissue_eqtl.pdf"},
        content = function(file) {
            pdf(file, width = 10, height = 8)
            
            power_unbalanced <- array(numeric(sizeLen.data()*mafLen1), dim=c(sizeLen.data(),mafLen1), dimnames = list(sizeVec.data(), as.character(mafVec1*100)))
            
            for (i in 1:sizeLen.data()){
                for (j in 1:mafLen1){
                    # treating each cell independantly 
                    power_unbalanced[i,j] <- powerEQTL(MAF=mafVec1[j], alpha=input$alpha1, myntotal=sizeVec.data()[i], delta = input$delta)
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
            mtext(paste("eQTL p-value = ", input$alpha1, " (one-way unbalanced ANOVA)", sep = ""),3)
            abline(v=0, h=seq(0,1,.1), lty=2, col="grey89")
            abline(h=0, v=c(1:10), lty=2,col="grey89")
            # add power curves
            for (i in 1:sizeLen.data()){
                lines(mafVec1*100, power_unbalanced[i,], type="l", lwd=4, col=colors[i])
                pos=which(power_unbalanced[,i]>=0.8)[1]
                abline(v=mafVec1[pos], col=colors[i])
                text(mafVec1[pos],0.8,labels = as.character(mafVec1[pos]), srt = 90,cex=0.6, adj=0.5)
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
                    helpText("$$\\text{= # of folds * standard deviation}$$"),
                    helpText("$$\\text{family-wise type I error rate = 0.05/nTests}$$"),
                    helpText("$$\\text{nTests = number of genes * number of SNPs}$$"))
    })
    
    output$description <- renderUI({
        withMathJax(
            helpText("This is more realistic model when treating neurons from one subject as an intra-class. "),
            helpText("We assume the following linear mixed effects model to characterize the association between genotype and gene expression:"),
            helpText(" $$y_{ij} = \\beta_{0i} + \\beta_1  x_i + \\epsilon_{ij}$$"),
            helpText(" $$\\beta_{0i} \\text{ ~ } N(\\beta_0, \\sigma^2_{\\beta})$$ "),
            helpText(" $$\\epsilon_{ij} \\text{ ~ } N(0, \\sigma^2_{\\epsilon})$$"),
            helpText("$$i = 1 ... n$$"),
            helpText("$$j = 1 ... m$$"))
    })
    
    output$description3 <- renderUI({
        withMathJax(
            helpText("$$\\text{family-wise type I error rate = 0.05/nTests}$$"),
            helpText("$$\\text{nTests = number of genes * number of SNPs}$$"))
    })
    
    output$description2 <- renderUI({
        withMathJax(
            helpText("This is a model when treating each tissue from one subject as an intra-class. "),
            helpText("We assume the following one-way unbalanced anova model based on Dupont and Plummer (1998):"),
            helpText(" $$y_{i} = \\gamma + \\lambda  x_i + \\epsilon_{i}$$"),
            helpText(" $$\\epsilon_{i} \\text{ ~ } N(0, \\rho^2)$$"),
            helpText("$$i = 1 ... n$$"))
    })
    
    
    ## Approximate sample sizes with given inputs
    
    ## First start from an initial guess: n = 100
    ## Then increment n by 100 in each iteration to approximate n when power level approaches input with error < 0.001
    ## Each time power level > input, cut increment value by 1/2
    ## Power level shall be < 1.0
    output$approx <- renderUI({
        # n & inc & approx for single cell eqtl
        # n2 & inc2 & approx2 for tissue eqtl
        n <- 100
        n2 <- 100
        inc <- 100
        inc2 <- 100
        approx <- power.eQTL.scRNAseq(delta=input$delta2, n=n, m=input$m2,sigma.y=input$sigma2, theta=input$maf2, rho=text1.data(), alpha=input$alpha2)
        approx2 <- powerEQTL(MAF=input$maf2, alpha=input$alpha2, myntotal=n2, delta = input$delta2)
        
        while (isTRUE(approx - input$power < 0 || approx - input$power > 0.005)) 
        {
            if (approx > input$power)
            {
                n = n - inc
                inc = inc/2
            }
            n = n + inc
            approx = power.eQTL.scRNAseq(delta=input$delta2, n=n, m=input$m2,sigma.y=input$sigma2, theta=input$maf2, rho=text1.data(), alpha=input$alpha2)
        }
        
        while (isTRUE(approx2 - input$power < 0 || approx2 - input$power > 0.005)) 
        {
            if (approx2 > input$power)
            {
                n2 = n2 - inc2
                inc2 = inc2/2
            }
            n2 = n2 + inc2
            approx2 = powerEQTL(MAF=input$maf2, alpha=input$alpha2, myntotal=n2, delta = input$delta2)
        }
        
        withMathJax(
            helpText("The minium sample size to approximate power level ", input$delta2, "(single-cell eQTL) is ", ceiling(n)),
            helpText("The minium sample size to approximate power level ", input$delta2, "(tissue eQTL) is ", ceiling(n2))
        )

    })
    
    
}