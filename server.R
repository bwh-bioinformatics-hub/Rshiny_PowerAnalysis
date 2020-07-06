library(shiny)
library(shinyjs)
library(devtools)

source("./global.R")

function(input, output, session) {
    
    ## Global Values
    
    ## values for theta (MAF) (single-cell)
    mafVec <- reactive({
        # read in maf range
        rg <- as.numeric(unlist(strsplit(as.character(input$MAF_sc), ",")))
        
        # check input values
        shiny::validate(need((length(rg) == 2), "Need one minimum and one maximum for MAF range")) 
        shiny::validate(need((rg[1] <= 0.5 && rg[1] >= 0 && rg[2] <= 0.5 && rg[2] >= 0), 
                             "Minor Allele Frequencies must in range [0, 0.5]"))
        shiny::validate(need((rg[1] <= rg[2]), "Max MAF value must be smaller than min MAF value"))
        
        # generate maf sequence (log-spaced)
        exp(seq(from=log(rg[1]), to=log(rg[2]), length.out = 100))
        })
    mafLen <-reactive(length(mafVec()))
    # values for theta (MAF) (tissue)
    mafVec1 <- reactive({
        # read in maf range
        rg1 <- as.numeric(unlist(strsplit(as.character(input$MAF_t), ",")))
        
        # check input values
        shiny::validate(need((length(rg1) == 2), "Need one minimum and one maximum for MAF range"))
        shiny::validate(need((rg1[1] <= 0.5 && rg1[1] >= 0 && rg1[2] <= 0.5 && rg1[2] >= 0), 
                             "Minor Allele Frequencies must in range [0, 0.5]"))
        shiny::validate(need((rg1[1] <= rg1[2]), "Max MAF value must be smaller than min MAF value"))
        
        # generate maf sequence (log-spaced)
        exp(seq(from=log(rg1[1]), to=log(rg1[2]), length.out = 100))
    })
    mafLen1 <-reactive(length(mafVec1()))
    
    ## values for sample sizes (tissue eqtl)
    # subVec.data <- reactive({unlist(strsplit(as.character(input$subjects), ","))})
    sizeVec.data <- reactive({as.numeric(unlist(strsplit(as.character(input$myntot), ",")))})
    sizeLen.data <- reactive(length(sizeVec.data()))
    # values for sample sizes (single-cell eqtl)
    # subVec1.data <- reactive({unlist(strsplit(as.character(input$subjects1), ","))})
    sizeVec1.data <- reactive({as.numeric(unlist(strsplit(as.character(input$myntot1), ",")))})
    sizeLen1.data <- reactive(length(sizeVec1.data()))
    
    ## validate intra-class correlation
    text.data <- reactive({
         shiny::validate(need((input$rho <= 1 && input$rho >= 0),
                       "Intra-class correlation must in range [0, 1]"))
         input$rho
    })    #sc
    text1.data <- reactive({
        shiny::validate(need((input$rho2 <= 1 && input$rho2 >= 0),
                             "Intra-class correlation must in range [0, 1]"))
        input$rho2
    })   #t
    
    ##initialize power
    power_sc = 0  ## (sc)
    power_t = 0 ##(t)
    
    ## sc power eQTL
    power_sc_eQTL <- reactive({
        # calculate scRNAseq for each sample size and maf
        sc = array(numeric(sizeLen1.data()*mafLen()), dim=c(sizeLen1.data(),mafLen()), 
                                    dimnames = list(as.character(sizeVec1.data()), as.character(mafVec())))
        for(i in 1:sizeLen1.data()) {
            for(j in 1:mafLen()) {
                sc[i,j]=power.eQTL.scRNAseq(delta=input$delta, n=sizeVec1.data()[i], 
                                                       m=input$m,sigma.y=input$sigma, theta=mafVec()[j], 
                                                       rho=text.data(), alpha=input$alpha/input$nTest)
            }
        }
        # store it to power_sc_eQTL
        sc
        })
    
    ## tissue power eQTL
    power_tissue <- reactive({
        t <- array(numeric(sizeLen.data()*mafLen1()), dim=c(sizeLen.data(),mafLen1()), 
                   dimnames = list(sizeVec.data(), as.character(mafVec1()*100)))
        
        # ANOVA 
        if (input$radio == "One-way unbalanced anova")
        {
            for (i in 1:sizeLen.data()){
                for (j in 1:mafLen1()){
                    # treating each cell independantly 
                    t[i,j] <- powerEQTL.ANOVA(MAF=mafVec1()[j], FWER=input$alpha1, nTests=input$nTest1, 
                                              n=sizeVec.data()[i], deltaVec = c(input$delta1, input$delta2), sigma = input$sigma1)
                }
            }
        } 
        else
        {
            # SLR
            for (i in 1:sizeLen.data()){
                for (j in 1:mafLen1()){
                    # treating each cell independantly 
                    t[i,j] <- powerEQTL.SLR.default(MAF=mafVec1()[j], FWER=input$alpha1, nTests = input$nTest1, n=sizeVec.data()[i], slope = input$delta1, sigma.y = input$sigma1)
                }
            }
        }
        
        # store it to power_tissue
        t
    })
    
    ##author names
    output$authors <- renderUI({
        tags$div(
            tags$h5("This R shiny application is developed by Xiaoqi Li, under the supervisions of Drs. 
                    Xianjun Dong and Weiliang Qiu."),
            
            paste0("Thank you very much for using our 'Sample size and power calculator for bulk tissue 
                   and single-cell eQTL analysis'! Please let us know what do you think."),
            tags$a("Click here!", href = "https://gitreports.com/issue/bwh-bioinformatics-hub/Rshiny_PowerAnalysis")
        )
    })
    
    
    ## single-cell eQTL
    
    ## hidden input options
    observeEvent(input$btn_sc, {
        toggle("MAF_sc")
        toggle("pos")
    })
    ## plot single-cell eQTL
    output$cell <- renderPlot({
        xrange <- range(mafVec()*100)
        yrange <- c(0:1)
        colors <- hcl.colors(sizeLen.data(),"Set 2")
        title <- "Single-cell eQTL"
        
        
        # treating each cell independantly 
        plot(xrange, yrange, xlim=range(xrange), log='x', type="n",
             xlab="MAF (%)",
             ylab="Power",
             main= title, xaxt = "n", yaxt = "n")
        axis(1, cex.axis = 1.2)
        axis(2, cex.axis = 1.2)
        
        
        mtext(paste("linear mixed effects model (alpha=",input$alpha/input$nTest,", slope=",input$delta,", n_cell="
                    , input$m ,", sigma=", input$sigma ,", rho=", input$rho ,")", sep = ""),3)
        
        abline(v=0, h=seq(0,1,.1), lty=2, col="grey89")
        abline(h=0, v=c(1:10), lty=2,col="grey89")

        # add power curves
        for (i in 1:sizeLen1.data()){
            lines(mafVec()*100, power_sc_eQTL()[i,], type="l", lwd=4, col=colors[i])
        }
        legend(input$pos, title="Sample size (n)", paste("n = ",sizeVec1.data(),"", sep=""),
               title.col='black',text.col=colors,cex =1, bty='n')
    })
    ## Interactive plot for sc
    observeEvent(input$sc_plot_hover, {
        # get hover coordinate
        power_sc <<- input$sc_plot_hover$y

        #redraw the sc power curves
        output$cell <- renderPlot({
            
            xrange <- range(mafVec()*100)
            yrange <- c(0:1)
            colors <- hcl.colors(sizeLen.data(),"Set 2")
            title <- "Single-cell eQTL"
            
            
            # treating each cell independantly 
            plot(xrange, yrange, xlim=range(xrange), log='x', type="n",
                 xlab="MAF (%)",
                 ylab="Power",
                 main= title, xaxt = "n", yaxt = "n")
            axis(1, cex.axis = 1.2)
            axis(2, cex.axis = 1.2)
            
            
            mtext(paste("linear mixed effects model (FWER=",input$alpha, ", nTests=", input$nTest,", slope=",input$delta,", n_cell="
                        , input$m ,", sigma=", input$sigma ,", rho=", input$rho ,")", sep = ""),3)
            
            abline(v=0, h=seq(0,1,.1), lty=2, col="grey89")
            abline(h=0, v=c(1:10), lty=2,col="grey89")
            abline(h=power_sc, lty=2, lwd = 1.5, col = 1)
            text(mafVec()[1]*100, power_sc,labels = as.character(round(power_sc,2)), cex=1, adj=0.5, pos = 1, offset = -.8)  #power level

            # add power curves
            for (i in 1:sizeLen1.data()){
                lines(mafVec()*100, power_sc_eQTL()[i,], type="l", lwd=4, col=colors[i])
                pos=which(power_sc_eQTL()[i,]>=power_sc)[1] #find the least maf for that power
                abline(v=mafVec()[pos]*100, col=colors[i]) #horizontal line for power level
                text(mafVec()[pos]*100,power_sc,labels = as.character(round(mafVec()[pos]*100,2)), srt = 90,cex=1, adj=0.5) #maf value
            }
            
            legend(input$pos, title="Sample size (n)", paste("n = ",sizeVec1.data(),"", sep=""),
                   title.col='black',text.col=colors,cex =1, bty='n')
        })
        
    })
    
    
    ## tissue eQTL
    
    ## hidden input options
    observeEvent(input$btn_t, {
        toggle("MAF_t")
        toggle("pos_t")
    })
    ## plot tissue eQTL
    output$tissue <- renderPlot({
        # ANOVA 
        if (input$radio == "One-way unbalanced anova")
        {
            ## plot
            
            # set up graph
            xrange <- range(mafVec1()*100)
            yrange <- c(0:1)
            colors <- hcl.colors(sizeLen.data(),"Set 2")
            title <- "Tissue eQTL"
            
            ## 	 # treating each cell independantly 
            plot(xrange, yrange, xlim=range(xrange), log='x', type="n",
                 xlab="MAF (%)",
                 ylab="Power",
                 main= title)
            mtext(paste("One-way unbalanced ANOVA (FWER=",input$alpha1, ", nTests=", input$nTest1,", Δ=",input$delta1,", sigma=", input$sigma1 ,")", sep = ""),3)
            
        } 
        # SLR
        else
        {
            ## plot
            
            # set up graph
            xrange <- range(mafVec1()*100)
            yrange <- c(0:1)
            colors <- hcl.colors(sizeLen.data(),"Set 2")
            title <- "Tissue eQTL"
            
            ## 	 # treating each cell independantly 
            plot(xrange, yrange, xlim=range(xrange), log='x', type="n",
                 xlab="MAF (%)",
                 ylab="Power",
                 main= title)
            mtext(paste("Simple Linear Regression (FWER=",input$alpha, ", nTests=", input$nTest1,", β=",input$delta1,", sigma=", input$sigma1 ,")", sep = ""),3)
        }
        
        abline(v=0, h=seq(0,1,.1), lty=2, col="grey89")
        abline(h=0, v=c(1:10), lty=2,col="grey89")
        
        # add power curves
        for (i in 1:sizeLen.data()){
            lines(mafVec1()*100, power_tissue()[i,], type="l", lwd=4, col=colors[i])
        }
        legend(input$pos_t, title="Sample size (n)", paste("n = ",sizeVec.data(), sep=""),
               title.col='black',text.col=colors,cex =1, bty='n')
    })
    # Reactive Input Title for delta for different models
    observeEvent(input$radio,{
        
        if(input$radio=="One-way unbalanced ANOVA")
        {
            updateNumericInput(session, inputId = "delta1", label  = "Mean difference of gene expression levels between genotype classes (δ)", value = input$delta1)
            show("delta2")
            output$demo <- renderUI({
                img(src='eQTL_demo_ANOVA.png', align = "middle", width = "550px", height = "170px")
            })
            output$description2 <- renderUI({
                tags$div(
                    tags$h5("One-way Unbalanced ANOVA Model"),
                    
                    withMathJax(),
                    
                    paste0("Power calculation for eQTL analysis that tests if a SNP is associated to a gene expression level by using unbalanced one-way ANOVA, as the 
                   GTEx consortium used (GTEx consortium, Nature Genetics, 2013). Assuming we are testing ", input$nTest1 , " SNP-gene eQTL pairs (i.e. number of 
                   genes * number of SNPs), using a Bonferroni correction, we set the significance threshold, α, to be ", input$alpha1, "/", input$nTest1, ". 
                   We model the expression data as log-normally distributed with a log standard deviation of ", input$sigma1, " (we assume that each genotype 
                   class of subjects have the same standard deviation). This level of noise can be based on estimates from previous study or pilot data. The effect size 
                   depends both on the minor allele frequency (MAF) of the SNP and the actual log expression change between genotype classes (denoted by δ).
                   Figure above shows the statistical power of an eQTL analysis using an ANOVA statistical test as a function the number of subjects and 
                   the minor allele frequency (MAF), and assumes δ=", input$delta1," (equivalent to detecting a log-expression change similar to the ", 
                           input$delta1/input$sigma1, "-fold of standard deviation within a single genotype class)."),
                    tags$br(),
                    paste0("According to"),
                    tags$a( "SAS online document", href =  "https://support.sas.com/documentation/cdl/en/statug/63033/HTML/default/viewer.htm#statug_power_a0000000982.htm"),
                    paste0(", the power calculation formula is:"),
                    helpText(" $$power = Pr(F \\ge F_{1-\\alpha}(k-1, N-k) | F\\text{ ~ }F_{k-1, N-k, \\lambda})$$"),
                    HTML(paste("where k = 3 is the number of groups of subjects, N is the total number of subjects. F",tags$sub("1−α"),"(k − 1, N − k) is the 100(1 − α)-th 
                   percentile of central F distribution with degrees of freedoms k − 1 and N − k, and F",tags$sub("k−1,N−k,λ")," is the non-central F distribution
                   with degrees of freedoms k − 1 and N − k and non-central parameter (ncp) λ. The ncp λ is equal to ")),
                    tags$br(),
                    helpText("$$\\lambda = N/\\sigma^2 (p^2q(1+p)\\delta_1^2 + q^2p(1+q)\\delta_2^2 + 2p^2q^2\\delta_1\\delta_2)$$"),
                    HTML(paste("when δ",tags$sub("1"),"equals to δ",tags$sub("2"),", we have")),
                    helpText(" $$\\lambda = 2 p (1-p)  (\\frac{\\delta}{\\sigma})$$"),
                    paste0("where p is the minor allele frequency, δ is the mean difference of gene expression levels between genotype classes, and σ is 
                   the standard deviation of gene expression levels in one group of subjects."),
                    tags$br(),
                    
                    paste0("This plot uses the function of ‘powerEQTL.ANOVA’ in R package "),
                    tags$a( "‘powerEQTL’", href = "https://CRAN.R-project.org/package=powerEQTL"),
                    paste0(". More details can be found in the online manual of the function ‘powerEQTL.ANOVA’ in R package "),
                    tags$a( "‘powerEQTL’", href = "https://CRAN.R-project.org/package=powerEQTL"),
                    tags$br(),
                    tags$br(),
                    
                    tags$b("References"),
                    tags$br(),
                    
                    paste0("The GTEx Consortium. The Genotype-Tissue Expression (GTEx) project. Nature Genetics. 2013 June ; 45(6): 580–585")
                    
                )
            })
        }
        else
        {
            updateNumericInput(session, inputId = "delta1", label = "Slope of the regression line (β)", value = input$delta1)
            hide("delta2")
            output$demo <- renderUI({
                img(src='eQTL_demo_SLR.png', align = "middle", width = "550px", height = "170px")
            })
            output$description2 <- renderUI({
                tags$div(
                    
                    withMathJax(),
                    
                    tags$h5("Simple Linear Regression"),

                    paste("To test if a SNP is associated with a gene expression level, we use the simple linear regression based on Dupont and Plummer (1998):"),
                    helpText(" $$y_{i} = \\beta_0 + \\beta_1  x_i + \\epsilon_{i}$$"),
                    helpText("$$i = 1 ... n$$"),
                    HTML(paste0("where y", tags$sub("i"), " is the expression level of the gene for the subject i and x", tags$sub("i"), " is the genotype of the i-th 
                    subject by using additive coding (e.g. 0 for AA, 1 for AC, and 2 for CC). β", tags$sub("1")," is the slope of the regression line and can be 
                    estimated with the following distribution under the alternative hypothesis:")),
                    tags$br(),
                    helpText(" $$\\beta_1 \\text{ ~ } N(\\delta, \\frac{\\rho^2}{L_{xx}}) \\text{ and } \\beta_0 = y - \\beta_1 x$$"),
                    paste0("Hence, the power can be estimated with: "),
                    helpText("$$1-\\beta = 1 - T_{n-2, \\lambda}[t_{n-2}(\\alpha/2)] + T_{n-2, \\lambda}[-t_{n-2}(\\alpha/2)]$$"),
                    HTML(paste0("where T", tags$sub("n−2,λ"), "(a) is the value at a of the cumulative distribution function of non-central t distribution with (n − 2)
                    degrees of freedom and non-centrality parameter λ. λ can be written as: ")),
                    helpText("$$ \\lambda = \\frac{\\delta}{\\sqrt{(\\sigma_y^2 - \\delta^2 2(1-p)p) / ((n-1)2(1-p)p)}} $$"),
                    
                    
                    paste0("This plot uses the function of ‘powerEQTL.SLR’ in R package "),
                    tags$a( "‘powerEQTL’", href = "https://CRAN.R-project.org/package=powerEQTL"),
                    paste0(". More details can be found in the online manual of the function ‘powerEQTL.SLR’ in R package "),
                    tags$a( "‘powerEQTL’", href = "https://CRAN.R-project.org/package=powerEQTL"),
                    tags$br(),
                    tags$br(),
                    
                    tags$b("References"),
                    tags$br(),
                    
                    paste0("Dupont, W.D. and Plummer, W.D.. Power and Sample Size Calculations for Studies Involving Linear Regression. Controlled Clinical Trials. 1998;19:589-601.")
                )
            })
        }
        ## update plot
        output$tissue <- renderPlot({
            power_tissue <- array(numeric(sizeLen.data()*mafLen1()), dim=c(sizeLen.data(),mafLen1()), dimnames = list(sizeVec.data(), as.character(mafVec1()*100)))
            
            # ANOVA 
            if (input$radio == "One-way unbalanced ANOVA")
            {
                for (i in 1:sizeLen.data())
                {
                    for (j in 1:mafLen1())
                    {
                        # treating each cell independantly 
                        power_tissue[i,j] <- powerEQTL.ANOVA(MAF=mafVec1()[j], FWER=input$alpha1, nTests=input$nTest1, 
                                                             n=sizeVec.data()[i], deltaVec = c(input$delta1, input$delta2), sigma = input$sigma1)
                    }
                }
                
                ## plot
                
                # set up graph
                xrange <- range(mafVec1()*100)
                yrange <- c(0:1)
                colors <- hcl.colors(sizeLen.data(),"Set 2")
                title <- "Tissue eQTL"
                
                ## 	 # treating each cell independantly 
                plot(xrange, yrange, xlim=range(xrange), log='x', type="n",
                     xlab="MAF (%)",
                     ylab="Power",
                     main= title)
                mtext(paste("One-way unbalanced ANOVA (FWER=",input$alpha1, ", nTests=", input$nTest1,", Δ=",input$delta1,", sigma=", input$sigma1 ,")", sep = ""),3)
                
            } 
            # SLR
            else
            {
                for (i in 1:sizeLen.data())
                {
                    for (j in 1:mafLen1())
                    {
                        power_tissue[i,j] <- powerEQTL.SLR.default(MAF=mafVec1()[j], FWER=input$alpha1, nTests = input$nTest1, n=sizeVec.data()[i], slope = input$delta1, sigma.y = input$sigma1)
                    }
                }

                # set up graph
                xrange <- range(mafVec1()*100)
                yrange <- c(0:1)
                colors <- hcl.colors(sizeLen.data(),"Set 2")
                title <- "Tissue eQTL"
                
                ## 	 # treating each cell independantly 
                plot(xrange, yrange, xlim=range(xrange), log='x', type="n",
                     xlab="MAF (%)",
                     ylab="Power",
                     main= title)
                mtext(paste("Simple Linear Regression (FWER=",input$alpha, ", nTests=", input$nTest1,", β=",input$delta1,", sigma=", input$sigma1 ,")", sep = ""),3)
            }
            
            abline(v=0, h=seq(0,1,.1), lty=2, col="grey89")
            abline(h=0, v=c(1:10), lty=2,col="grey89")
            
            # add power curves
            for (i in 1:sizeLen.data()){
                lines(mafVec1()*100, power_tissue[i,], type="l", lwd=4, col=colors[i])
            }
            legend(input$pos_t, title="Sample size (n)", paste("n = ",sizeVec.data(), sep=""),
                   title.col='black',text.col=colors,cex =1, bty='n')
        })
    })
    ## Interactive plot for tissue
    observeEvent(input$t_plot_hover$y,{
        # get hover coordinate
        power_t <<- input$t_plot_hover$y
        
        # ANOVA 
        if (input$radio == "One-way unbalanced ANOVA")
        {
            ## replot tissue eQTL for ANOVA
            output$tissue <- renderPlot({
                # set up graph
                xrange <- range(mafVec1()*100)
                yrange <- c(0:1)
                colors <- hcl.colors(sizeLen.data(),"Set 2")
                title <- "Tissue eQTL"
                
                ## 	 # treating each cell independantly 
                plot(xrange, yrange, xlim=range(xrange), log='x', type="n",
                     xlab="MAF (%)",
                     ylab="Power",
                     main= title)
                mtext(paste("One-way unbalanced ANOVA (FWER=",input$alpha1, ", nTests=", input$nTest1,", Δ=",input$delta1,", sigma=", input$sigma1 ,")", sep = ""),3)
                
                abline(v=0, h=seq(0,1,.1), lty=2, col="grey89")
                abline(h=0, v=c(1:10), lty=2,col="grey89")
                abline(h=power_t, lty=2, lwd = 1.5, col = 1)
                text(mafVec1()[1]*100, power_t,labels = as.character(round(power_t,2)), cex=1, adj=0.5, pos = 1, offset = -.8)  #power level
                
                # add power curves
                for (i in 1:sizeLen.data()){
                    lines(mafVec1()*100, power_tissue()[i,], type="l", lwd=4, col=colors[i])
                    pos1=which(power_tissue()[i,]>=power_t)[1]
                    abline(v=mafVec1()[pos1]*100, col=colors[i])
                    text(mafVec1()[pos1]*100,power_t,labels = as.character(round(mafVec1()[pos1]*100,2)), srt = 90,cex=1, adj=0.5)
                }
                
                legend(input$pos_t, title="Sample size (n)", paste("n = ",sizeVec.data(), sep=""),
                       title.col='black',text.col=colors,cex =1, bty='n')
            })
        } 
        # SLR
        if (input$radio == "Simple linear regression")
        {
            ## plot
            output$tissue <- renderPlot({
                # set up graph
                xrange <- range(mafVec1()*100)
                yrange <- c(0:1)
                colors <- hcl.colors(sizeLen.data(),"Set 2")
                title <- "Tissue eQTL"
                
                ## 	 # treating each cell independantly 
                plot(xrange, yrange, xlim=range(xrange), log='x', type="n",
                     xlab="MAF (%)",
                     ylab="Power",
                     main= title)
                mtext(paste("Simple Linear Regression (FWER=",input$alpha1, ", nTests=", input$nTest1,", β=",input$delta1,", sigma=", input$sigma1 ,")", sep = ""),3)
                
                abline(v=0, h=seq(0,1,.1), lty=2, col="grey89")
                abline(h=0, v=c(1:10), lty=2,col="grey89")
                abline(h=power_t, lty=2, lwd = 1.5, col = 1)
                text(mafVec1()[1]*100, power_t,labels = as.character(round(power_t,2)), cex=1, adj=0.5, pos = 1, offset = -.8)  #power level
                
                # add power curves
                for (i in 1:sizeLen.data()){
                    lines(mafVec1()*100, power_tissue()[i,], type="l", lwd=4, col=colors[i])
                    pos1=which(power_tissue()[i,]>=power_t)[1]
                    abline(v=mafVec1()[pos1]*100, col=colors[i])
                    text(mafVec1()[pos1]*100,power_t,labels = as.character(round(mafVec1()[pos1]*100,2)), srt = 90,cex=1, adj=0.5)
                }
                legend(input$pos_t, title="Sample size (n)", paste("n = ",sizeVec.data(), sep=""),
                       title.col='black',text.col=colors,cex =1, bty='n')
            })
        }
            
        })


    ## clicking on the export button will generate a pdf file for sc eQTL
    output$export = downloadHandler(
        filename = function() {"single_cell_eqtl.pdf"},
        content = function(file) {
            pdf(file, width = 8, height = 6)

            xrange <- range(mafVec()*100)
            yrange <- c(0:1)
            colors <- hcl.colors(sizeLen.data(),"Set 2")
            title <- "Single-cell eQTL"
            
            
            # treating each cell independantly 
            plot(xrange, yrange, xlim=range(xrange), log='x', type="n",
                 xlab="MAF (%)",
                 ylab="Power",
                 main= title)
            
            mtext(paste("linear mixed effects model (FWER=",input$alpha, ", nTests=", input$nTest,", slope=",input$delta,", n_cell="
                        , input$m ,", sigma=", input$sigma ,", rho=", input$rho ,")", sep = ""),3)
            
            abline(v=0, h=seq(0,1,.1), lty=2, col="grey89")
            abline(h=0, v=c(1:10), lty=2,col="grey89")
            if (power_sc != 0)
            {
                abline(h=power_sc, lty=2, lwd = 1.5, col = 1)
                text(mafVec()[1]*100, power_sc,labels = as.character(round(power_sc,2)), cex=1, adj=0.5, pos = 1, offset = -.8)  #power level
            }
            
            # add power curves
            for (i in 1:sizeLen1.data())
            {
                lines(mafVec()*100, power_sc_eQTL()[i,], type="l", lwd=4, col=colors[i])
                if (power_sc != 0)
                {
                    pos=which(power_sc_eQTL()[i,]>=power_sc)[1]
                    abline(v=mafVec()[pos]*100, col=colors[i])
                    text(mafVec()[pos]*100,power_sc,labels = as.character(round(mafVec()[pos]*100,2)), srt = 90,cex=1, adj=0.5)
                }
            }
            
            legend(input$pos, title="Sample size (n)", paste("n = ",sizeVec.data(), sep=""),
                   title.col='black',text.col=colors,cex =1, bty='n')
            
            dev.off()
        }
    )
    ## clicking on the export button will generate a pdf file for tissue eQTL
    output$export2 = downloadHandler(
        filename = function() {"tissue_eqtl.pdf"},
        content = function(file) {
            pdf(file, width = 8, height = 6)
            
            # ANOVA 
            if (input$radio == "One-way unbalanced ANOVA")
            {
                ## plot
                
                # set up graph
                xrange <- range(mafVec1()*100)
                yrange <- c(0:1)
                colors <- hcl.colors(sizeLen.data(),"Set 2")
                title <- "Tissue eQTL"
                
                ## 	 # treating each cell independantly 
                plot(xrange, yrange, xlim=range(xrange), log='x', type="n",
                     xlab="MAF (%)",
                     ylab="Power",
                     main= title)
                mtext(paste("One-way unbalanced ANOVA (FWER=",input$alpha1, ", nTests=", input$nTest1,", Δ=",input$delta1,", sigma=", input$sigma1 ,")", sep = ""),3)
                
            } 
            # SLR
            else
            {
                ## plot
                
                # set up graph
                xrange <- range(mafVec1()*100)
                yrange <- c(0:1)
                colors <- hcl.colors(sizeLen.data(),"Set 2")
                title <- "Tissue eQTL"
                
                ## 	 # treating each cell independantly 
                plot(xrange, yrange, xlim=range(xrange), log='x', type="n",
                     xlab="MAF (%)",
                     ylab="Power",
                     main= title)
                mtext(paste("Simple Linear Regression (FWER=",input$alpha1, ", nTests=", input$nTest1,", β=",input$delta1,", sigma=", input$sigma1 ,")", sep = ""),3)
            }
            
            abline(v=0, h=seq(0,1,.1), lty=2, col="grey89")
            abline(h=0, v=c(1:10), lty=2,col="grey89")
            if (power_t != 0)
            {
                abline(h=power_t, lty=2, lwd = 1.5, col = 1)
                text(mafVec1()[1]*100, power_t,labels = as.character(round(power_t,2)), cex=1, adj=0.5, pos = 1, offset = -.8)  #power level
            }
                        
            # add power curves
            for (i in 1:sizeLen.data())
            {
                lines(mafVec1()*100, power_tissue()[i,], type="l", lwd=4, col=colors[i])
                if (power_t != 0)
                {
                    pos1=which(power_tissue()[i,]>=power_t)[1]
                    abline(v=mafVec1()[pos1]*100, col=colors[i])
                    text(mafVec1()[pos1]*100,power_t,labels = as.character(round(mafVec1()[pos1]*100,2)), srt = 90,cex=1, adj=0.5)
                }
            }
            legend(input$pos_t, title="Sample size (n)", paste("n = ",sizeVec.data(), sep=""),
                   title.col='black',text.col=colors,cex =1, bty='n')
            
            dev.off()
        }
    )
    
    ## Sc Plot Explanations
    output$description <- renderUI({
        tags$div(
        
        withMathJax(),
    
        paste0("We assume the following simple linear mixed effects model for each (SNP, gene) pair to 
               characterize the association between genotype and gene expression:"),
        helpText(" $$y_{ij} = \\beta_{0i} + \\beta_1  x_i + \\epsilon_{ij}$$"),
        helpText(" $$\\beta_{0i} \\text{ ~ } N(\\beta_0, \\sigma^2_{\\beta})$$ "),
        helpText(" $$\\epsilon_{ij} \\text{ ~ } N(0, \\sigma^2_{\\epsilon})$$"),
        helpText("$$i = 1 ... n$$"),
        helpText("$$j = 1 ... m$$"),
        HTML(paste0("n is the number of subjects, m is the number of cells per subjectct, y", tags$sub("ij")," is 
        the gene expression level for the j-th cell of the i-th subject, x", tags$sub("i"),"is the genotype for 
        the i-th subject using additive coding. That is, x", tags$sub("i")," = 0 indicates the i-th subject is
        a wildtype homozygote, x", tags$sub("i")," = 1 indicates the i-th subject is a heterozygote, and x", 
        tags$sub("i")," = 2 indicates the i-th subject is a mutation homozygote.")),
        paste0("For a given SNP, we assume Hardy-Weinberg Equilibrium and denote the minor allele frequency of the SNP
        as θ. We can derive the power calculation formula is:"),
        helpText("$$power = 1 - \\Phi(z_{\\alpha*/2} - ab) + \\Phi(-z_{\\alpha*/2} - ab)$$"),
        helpText("$$\\text{where } a=\\frac{\\sqrt{2\\theta(1-\\theta)}}{\\beta} \\text{ and } b = \\frac{\\delta \\sqrt{m(n-1)}}{\\sqrt{1+(m-1)\\rho}}$$"),
        HTML(paste0("z",tags$sub("α∗/2"), " is the upper 100α∗/2 percentile of the standard normal distribution,
        α∗ = α/nTests, nTests is the number of (SNP, gene) pairs, ρ is the intraclass correlation,.")),
        
        tags$br(),
        paste0("This plot uses function ‘powerEQTL.scRNAseq’ in R package "),
        tags$a( "‘powerEQTL’", href = "https://CRAN.R-project.org/package=powerEQTL"),
        paste0(". More details can be found in the online manual of the function ‘powerEQTL.scRNAseq’ in R package "),
        tags$a( "‘powerEQTL’", href = "https://CRAN.R-project.org/package=powerEQTL"),
        tags$br(),
        tags$br(),   
        
        tags$b("References"),
        tags$br(),
        
        paste0("Dong X, Li X, Chang T, Weiss S, and Qiu W. powerEQTL: an R package and R shiny application for calculating 
               sample size and power of bulk tissue and single-cell eQTL analysis. manuscript. (2020)")
        )
    })
    
    ## Approximate sample sizes with given inputs
    # observeEvent(input$MAF_sc,{
    #     updateSliderInput()
    # })
    sc_est <- reactive({
        n=uniroot(f=diffPower4ss.ANOVA,
                  interval = c(10, 1e30),
                  MAF=input$maf_est,
                  deltaVec=c(input$slope_est, input$slope1_est),
                  power=input$power_est,
                  sigma=input$sigma_est,
                  FWER=input$FWER_est,
                  nTests=input$nTest_est
        )

        n2=uniroot(f=diffPower4ss.SLR,
                   interval = c(10, 1e30),
                   MAF=input$maf_est,
                   slope=input$slope_est,
                   power=input$power_est,
                   sigma=input$sigma_est,
                   FWER=input$FWER_est,
                   nTests=input$nTest_est
        )
        
        # updateSliderInput(inputId = "n_est", label = "Number_of_subjects_needed",
        #                   value = 100, min = 10, max = 1e30)
        data.frame(Model_used = c("One-way unbalanced Anova", "Simple Linear Regression"),
                   Number_of_subjects_needed = c(ceiling(n$root), ceiling(n2$root)),
                   Minor_Allele_Frequency = rep(input$maf_est,2),
                   Power_level = rep(input$power_est,2),
                   δ = c(paste0(input$slope_est, ", ", input$slope1_est),input$slope_est),
                   σ = rep(input$sigma_est,2)
        )
    })
    t_est <- reactive({
        n=uniroot(f=diffPower4ss.scRNAseq,
                  interval = c(1, 1e30),
                  MAF=input$maf_est,
                  m=input$m_est,
                  rho=input$rho_est,
                  slope=input$slope_est,
                  power=input$power_est,
                  sigma=input$sigma_est,
                  FWER=input$FWER_est,
                  nTests=input$nTest_est
        )
        
        data.frame(Model_used = c("Simple linear mixed effects model"),
                   No_subjects_needed = ceiling(n$root),
                   Minor_Allele_Frequency = input$maf_est,
                   Power_level = input$power_est,
                   slope = input$slope_est,
                   σ = input$sigma_est,
                   m = input$m_est,
                   ρ = input$rho_est
        )
    })
    output$approx <- renderTable(sc_est())
    output$approx1 <- renderTable(t_est())
    
    ##Explanation for Sample Estimation
    output$title <- renderUI(tags$h5("Tissue eQTL"))
    output$title1 <- renderUI(tags$div(tags$br(), tags$h5("Single-cell eQTL")))
    output$Explanation <- renderUI({
        tags$div(
            tags$br(),
            paste0("alpha = Family Type-I Error Rate / Number of Tests =", input$FWER_est," /", input$nTest_est," = ", input$FWER_est/input$nTest_est),
            tags$br(),
            paste0("eQTL effect size = mean difference of gene expression levels between groups(delta) / standard deviation of gene expression levels(sigma)"),
            tags$br(),
            tags$br(),
            paste0("These tables uses function ‘ssEQTL.ANOVA’, 'ssEQTL.scRNAseq', and 'ssEQTL.SLR' in R package "),
            tags$a( "‘powerEQTL’", href = "https://CRAN.R-project.org/package=powerEQTL"),
            paste0(". More details can be found in the online manual of the function ‘ssEQTL.ANOVA’, 'ssEQTL.scRNAseq', and 'ssEQTL.SLR' in R package "),
            tags$a( "‘powerEQTL’", href = "https://CRAN.R-project.org/package=powerEQTL"),
            tags$br(),
            tags$br(),
            
            tags$b("References"),
            tags$br(),
            
            paste0("Dong X, Li X, Chang T, Weiss S, and Qiu W. powerEQTL: an R package and R shiny application for calculating 
               sample size and power of bulk tissue and single-cell eQTL analysis. manuscript. (2020)"),
            tags$br(),
            paste0("Dupont, W.D. and Plummer, W.D.. Power and Sample Size Calculations for Studies Involving Linear Regression. Controlled Clinical Trials. 1998;19:589-601."),
            tags$br(),
            paste0("Lonsdale J and Thomas J, et al. The Genotype-Tissue Expression (GTEx) project. Nature Genetics, 45:580-585, 2013.")
            
            )
    })
}