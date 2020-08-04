library(shiny)
library(shinyjs)
library(devtools)

source("./global.R")
source("./package/ssEQTL.ANOVA.R")
source("./package/ssEQTL.SLR.R")
source("./package/ssEQTL.scRNAseq.R")

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
    mafLen <- reactive(length(mafVec()))
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
    mafLen1 <- reactive(length(mafVec1()))
    
    ## values for sample sizes (tissue eqtl)
    sizeVec.data <- reactive({as.numeric(unlist(strsplit(as.character(input$myntot1), ",")))})
    sizeLen.data <- reactive(length(sizeVec.data()))
    # values for sample sizes (single-cell eqtl)
    sizeVec1.data <- reactive({as.numeric(unlist(strsplit(as.character(input$myntot), ",")))})
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
    
    ## initialize power
    power_sc = 0.8  ## (sc)
    power_t = 0.8 ##(t)
    
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
    
    ## author names
    output$about <- renderUI({
        tags$div(
            withMathJax(),
            
            tags$br(),
            
            tags$h4("What is eQTL?"),
            tags$div(img(src='eQTL_demo.svg', align = "center", width = "550px", height = "170px"), style="text-align: center;"),
            tags$br(),
            
            paste0("eQTL stands for Expression Quantitative Trait Locus, which is a locus that explains a fraction of 
                   the genetic variation of a gene expression phenotype. Standard eQTL analysis involves a direct 
                   association test between markers of genetic variation with gene expression levels typically measured
                   in tens or hundreds of individuals."),
            tags$br(),
            tags$br(),
            
            paste0("Single-cell eQTL measures single-cell level gene expression across subjects. Tissue eQTL measures
                   tissue level gene expression across subjects."),
            tags$br(),
            tags$br(),
            
            
            tags$h4("One-way unbalanced ANOVA model"),

            paste0("Power calculation for eQTL analysis that tests if a SNP is associated to a gene expression level 
                    by using unbalanced one-way ANOVA, as the GTEx consortium used (GTEx consortium, Nature Genetics, 2013)."),
            tags$br(),
            tags$br(),
            HTML(paste0("Suppose there are k = 3 groups of subjects: (1) mutation homozygotes; (2) heterozygotes; and (3)
                        wildtype homozygotes. We would like to test if the mean expression &mu;<sub>i</sub>, i = 1, ..., k, of the
                        gene is the same among the k groups of subjects. We can use the following one-way ANOVA model to characterize
                        the relationship between observed gene expression level y<sub>ij</sub> and the population mean expression level
                        &mu;<sub>i</sub>:")),
            tags$br(),
            helpText(" $$y_{ij} = \\mu_i + \\epsilon_{ij}$$"),
            helpText(" $$\\epsilon_{ij} \\text{ ~ } N(0, \\sigma^2)$$"),
            helpText("$$i = 1 ... k$$"),
            helpText("$$j = 1 ... n$$"),
            HTML(paste0("where y<sub>ij</sub> is the observed gene expression level for the j-th subject in the i-th group, &mu;<sub>i</sub>
            is the mean gene expression level of the i-th group, &epsilon;<sub>ij</sub> is the random error, k is the number of groups, 
            n<sub>i</sub> is the number of subjects in the i-th group. Denote the total number of subjects as N = &Sigma;<sup>k</sup><sub>i=1</sub>
            n<sub>i</sub>. That is, we have n<sub>1</sub> mutation homozygotes, n<sub>2</sub> heterozygotes, and n<sub>3</sub> wildtype homozygotes.")),
            tags$br(),
            tags$br(),
            paste0("According to O'Brien and Muller (1993), the power calculation formula is:"),
            helpText(" $$power = Pr(F \\ge F_{1-\\alpha}(k-1, N-k) | F\\text{ ~ }F_{k-1, N-k, \\lambda})$$"),
            HTML(paste("where k = 3 is the number of groups of subjects, N is the total number of subjects. F",tags$sub("1−α"),"(k − 1, N − k) is the 100(1 − α)-th 
                   percentile of central F distribution with degrees of freedoms k − 1 and N − k, and F",tags$sub("k−1,N−k,λ")," is the non-central F distribution
                   with degrees of freedoms k − 1 and N − k and non-central parameter (ncp) λ. The ncp λ is equal to ")),
            tags$br(),
            helpText("$$\\lambda = \\frac{N}{\\sigma^2} (p^2q(1+p)\\delta_1^2 + q^2p(1+q)\\delta_2^2 + 2p^2q^2\\delta_1\\delta_2)$$"),
            HTML(paste("where δ",tags$sub("1")," represents the mean difference of gene expression between AA and AC (A is major allele and C is minor allele)
                        , and δ",tags$sub("2")," represents the mean difference of gene expression between CC and AC, p is the minor allele frequency, q = 1-p, 
                        and N is the sample size. ")),
            tags$br(),
            tags$br(),
            HTML(paste("When δ",tags$sub("1")," = δ",tags$sub("2"), " = δ, we have")),
            helpText(" $$\\lambda = 2 p q  N (\\frac{\\delta}{\\sigma})^2$$"),
            paste0("where p is the minor allele frequency, δ is the mean difference of gene expression levels between genotype classes, and σ is 
                   the standard deviation of random error in one group of subjects."),
            tags$br(),
            
            tags$br(),

            tags$h4("Simple linear regression"),

            paste("To test if a SNP is linearly associated with a gene expression level, we use the simple linear regression:"),
            helpText(" $$y_{i} = \\beta_0 + \\beta_1  x_i + \\epsilon_{i}$$"),
            helpText(" $$\\epsilon_{i} \\text{ ~ } N(0, \\sigma^2)$$"),
            helpText("$$i = 1 ... n$$"),
            HTML(paste0("where y", tags$sub("i"), " is the expression level of the gene for the subject i and x", tags$sub("i"), " is
                    the genotype of the i-th subject by using additive coding (e.g. 0 for AA, 1 for AC, and 2 for CC; A is major 
                    allele and C is minor allele), σ is the standard deviation of random error, and σ<sub>y</sub>
                    is the standard deviation of the outcome y<sub>i</sub>. ")),
            tags$br(),
            tags$br(),
            paste0("We can derive the power calculation formula as: "),
            helpText("$$power = 1 - T_{n-2, \\lambda}[t_{n-2}(\\alpha/2)] + T_{n-2, \\lambda}[-t_{n-2}(\\alpha/2)]$$"),
            HTML(paste0("where T", tags$sub("n−2,λ"), "(a) is the value at a of the cumulative distribution function of non-central t distribution with (n − 2)
                    degrees of freedom and non-centrality parameter λ. λ can be written as: ")),
            helpText("$$ \\lambda = \\frac{\\beta_1}{\\sqrt{(\\sigma_y^2 - \\beta_1^2 2pq) / [(n-1)2pq]}} $$"),
            paste0("p is the minor allele frequency (MAF) and q = 1-p."),
            
            tags$br(),
            
            tags$h4("Simple linear mixed effects model"),
            
            paste0("For single-cell eQTL, we assume the following simple linear mixed effects model for each (SNP, gene) pair to 
                    characterize the association between genotype and gene expression:"),
            helpText(" $$y_{ij} = \\beta_{0i} + \\beta_1  x_i + \\epsilon_{ij}$$"),
            helpText(" $$\\beta_{0i} \\text{ ~ } N(\\beta_0, \\sigma^2_{\\beta})$$ "),
            helpText(" $$\\epsilon_{ij} \\text{ ~ } N(0, \\sigma^2)$$"),
            helpText("$$i = 1 ... n$$"),
            helpText("$$j = 1 ... m$$"),
            HTML(paste0("n is the number of subjects, m is the number of cells per subjectct, y", tags$sub("ij")," is 
                        the gene expression level for the j-th cell of the i-th subject, x", tags$sub("i")," is the genotype for 
                        the i-th subject using additive coding. That is, x", tags$sub("i")," = 0 indicates the i-th subject is
                        a wildtype homozygote, x", tags$sub("i")," = 1 indicates the i-th subject is a heterozygote, and x", 
                        tags$sub("i")," = 2 indicates the i-th subject is a mutation homozygote, and σ<sub>y</sub>
                        is the standard deviation of the outcome y<sub>ij</sub>.")),
            tags$br(),
            tags$br(),
            paste0("For a given SNP, we assume Hardy-Weinberg Equilibrium and denote the minor allele frequency of the SNP
                as p. The power calculation formula is:"),
            helpText("$$power = 1 - \\Phi(z_{\\alpha*/2} - ab) + \\Phi(-z_{\\alpha*/2} - ab)$$"),
            helpText("$$\\text{where } a=\\frac{\\sqrt{2p(1-p)}}{\\sigma_y} \\text{ and } b = \\frac{\\beta_1 \\sqrt{m(n-1)}}{\\sqrt{1+(m-1)\\rho}}$$"),
            HTML(paste0("z",tags$sub("α∗/2"), " is the upper 100α∗/2 percentile of the standard normal distribution,
                        α∗ = α/nTests, nTests is the number of (SNP, gene) pairs, ρ is the intraclass correlation,.")),
            
            tags$br(),
            tags$br(),
            
            tags$b("References"),
            tags$br(),
            
            paste0("O'Brien, R. G., & Muller, K. E. (1993). Unified power analysis for t-tests through multivariate hypotheses.
                   In L. K. Edwards (Ed.), Statistics: Textbooks and monographs, Vol. 137. Applied analysis of variance in 
                   behavioral science (p. 297–344). Marcel Dekker."),
            tags$br(),
            paste0("The GTEx Consortium. The Genotype-Tissue Expression (GTEx) project. Nature Genetics, 45:580-585, 2013."),
            tags$br(),
            tags$br(),
            
            tags$h4("Developer and Bug Report"),
            paste0("This R shiny application is developed by Xiaoqi Li (xli778@wisc.edu), under the supervisions of Drs. 
                    Xianjun Dong and Weiliang Qiu."),
            tags$br(),
            
            paste0("Thank you very much for using our power calculators for bulk tissue 
                   and single-cell eQTL analysis! Please let us know what do you think by "),
            tags$a("clicking here!", href = "https://gitreports.com/issue/bwh-bioinformatics-hub/Rshiny_PowerAnalysis", target="_blank"),
            tags$br(),
            tags$br(),
            
            tags$h4("To cite our work, please use citation below:"),
            paste0("Dong X, Li X, Chang T, Weiss S, and Qiu W. powerEQTL: an R package and R shiny application for calculating 
                    sample size and power of bulk tissue and single-cell eQTL analysis. manuscript. (2020)"),
            tags$br(),
            tags$br(),
            
            tags$h4("Disclaimer"),
            tags$a(href="https://bioinformatics.bwh.harvard.edu/", "Genomics and Bioinformatics Hub ", target='_blank'),
            paste0(" cannot and will not be held legally, financially, or medically responsible for decisions made using its calculators, equations, content, and algorithms."),
            tags$br(),
            tags$br(),
            
        )
    })
    
    
    ## single-cell eQTL
    
    ## hidden input options
    observeEvent(input$btn_sc, {
        toggle("MAF_sc")
        toggle("pos")
        if (input$btn_sc %% 2 == 1) {
            txt <- "Less Options"
        } else {
            txt <- "More Options"
        }
        updateActionButton(session, "btn_sc", txt)
    })
    ## plot single-cell eQTL
    output$cell <- renderPlot({
        xrange <- range(mafVec()*100)
        yrange <- c(0:1)
        colors <- hcl.colors(sizeLen1.data(),"Set 2")
        par(mar = c(5, 5, 1, 2), cex.lab = 1.5)
        
        # treating each cell independantly 
        plot(xrange, yrange, xlim=range(xrange), log='x', type="n",
             xlab="MAF (%)",
             ylab="Power",
             xaxt = "n", yaxt = "n")
        axis(1, cex.axis = 1.5)
        axis(2, cex.axis = 1.5)

        abline(v=0, h=seq(0,1,.1), lty=2, col="grey89")
        abline(h=0, v=c(1:10), lty=2,col="grey89")
        abline(h=power_sc, lty=2, lwd = 1.5, col = 1)
        text(mafVec()[1]*100, power_sc,labels = as.character(round(power_sc,2)), cex=1.5, adj=0.5, pos = 1, offset = -1)  #power level
        
        # add power curves
        for (i in 1:sizeLen1.data()){
            lines(mafVec()*100, power_sc_eQTL()[i,], type="l", lwd=4, col=colors[i])
            pos=which(power_sc_eQTL()[i,]>=power_sc)[1] #find the least maf for that power
            abline(v=mafVec()[pos]*100, col=colors[i]) #horizontal line for power level
            text(mafVec()[pos]*100,power_sc,labels = as.character(round(mafVec()[pos]*100,2)), srt = 90,cex=1.5, adj=0.5) #maf value
        }
        
        legend(input$pos, title="Sample size (n)", legend = sizeVec1.data(),
               title.col='black', text.col=colors, cex =1.5, bty='n', pch = 15, col = colors)
    })
    ## Interactive plot for sc
    observeEvent(input$sc_plot_hover, {
        # get hover coordinate
        power_sc <<- input$sc_plot_hover$y

        #redraw the sc power curves
        output$cell <- renderPlot({
            
            xrange <- range(mafVec()*100)
            yrange <- c(0:1)
            colors <- hcl.colors(sizeLen1.data(),"Set 2")
            par(mar = c(5, 5, 1, 2), cex.lab = 1.5)
            
            # treating each cell independantly 
            plot(xrange, yrange, xlim=range(xrange), log='x', type="n",
                 xlab="MAF (%)",
                 ylab="Power",
                 xaxt = "n", yaxt = "n")
            axis(1, cex.axis = 1.5)
            axis(2, cex.axis = 1.5)
            
            abline(v=0, h=seq(0,1,.1), lty=2, col="grey89")
            abline(h=0, v=c(1:10), lty=2,col="grey89")
            abline(h=power_sc, lty=2, lwd = 1.5, col = 1)
            text(mafVec()[1]*100, power_sc, labels = as.character(round(power_sc,2)), cex=1.5, adj=0.5, pos = 1, offset = -1)  #power level

            # add power curves
            for (i in 1:sizeLen1.data()){
                lines(mafVec()*100, power_sc_eQTL()[i,], type="l", lwd=4, col=colors[i])
                pos=which(power_sc_eQTL()[i,]>=power_sc)[1] #find the least maf for that power
                abline(v=mafVec()[pos]*100, col=colors[i]) #horizontal line for power level
                text(mafVec()[pos]*100,power_sc,labels = as.character(round(mafVec()[pos]*100,2)), srt = 90,cex=1.5, adj=0.5) #maf value
            }
            
            legend(input$pos, title="Sample size (n)", legend = sizeVec1.data(),
                   title.col='black',text.col=colors,cex =1.5, bty='n', pch = 15, col = colors)
        })
        
    })
    
    
    ## tissue eQTL
    tissue <- reactive({
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
            
        }
        # set up graph
        xrange <- range(mafVec1()*100)
        yrange <- c(0:1)
        colors <- hcl.colors(sizeLen.data(),"Set 2")
        par(mar = c(5, 5, 1, 2), cex.lab = 1.5)
        
        ## 	 # treating each cell independantly 
        plot(xrange, yrange, xlim=range(xrange), log='x', type="n",
             xlab="MAF (%)",
             ylab="Power",
             xaxt = "n", yaxt = "n")
        
        axis(1, cex.axis = 1.5)
        axis(2, cex.axis = 1.5)
        
        abline(v=0, h=seq(0,1,.1), lty=2, col="grey89")
        abline(h=0, v=c(1:10), lty=2,col="grey89")
        abline(h=power_t, lty=2, lwd = 1.5, col = 1)
        text(mafVec1()[1]*100, power_t,labels = as.character(round(power_t,2)), cex=1.5, adj=0.5, pos = 1, offset = -1)  #power level
        
        # add power curves
        for (i in 1:sizeLen.data()){
            lines(mafVec1()*100, power_tissue[i,], type="l", lwd=4, col=colors[i])
            pos1=which(power_tissue[i,]>=power_t)[1]
            abline(v=mafVec1()[pos1]*100, col=colors[i])
            text(mafVec1()[pos1]*100,power_t,labels = as.character(round(mafVec1()[pos1]*100,2)), srt = 90,cex=1.5, adj=0.5)
        }
        
        legend(input$pos_t, title="Sample size (n)", legend = sizeVec.data(),
               title.col='black',text.col=colors,cex =1.5, bty='n', pch = 15, col = colors)
    })
    ## hidden input options
    observeEvent(input$btn_t, {
        toggle("MAF_t")
        toggle("pos_t")
        if (input$btn_t %% 2 == 1) {
            txt <- "Less Options"
        } else {
            txt <- "More Options"
        }
        updateActionButton(session, "btn_t", txt)
    })
    ## plot tissue eQTL
    # Reactive Input Title for delta for different models
    observeEvent(input$radio,{
        
        if(input$radio=="One-way unbalanced ANOVA")
        {
            updateNumericInput(session, inputId = "delta1", label  = "Mean difference of gene expression (δ1)", value = input$delta1)
            updateNumericInput(session, inputId = "sigma1", 
                               label = "Standard deviation of random error (σ)",
                               value = input$sigma1, min = 0.001, max = 1, step = 0.001)
            show("delta2")
            output$description2 <- renderUI({
                tags$div(
                    withMathJax(),
                    
                    # tags$br(),
                    
                    paste0("The figure above shows the statistical power of an eQTL analysis as a function of the minor 
                            allele frequency (MAF) with different sample sizes. Unbalanced one-way ANOVA test was used 
                            to test the potential nonlinear association of genotype classes to expression level, as Ref.
                            (GTEx consortium, Nature Genetics, 2013). We assume the random errors of preprocessed (e.g.,
                            log transformed) expression levels follow normal distribution with a standard deviation ",
                            input$sigma1, " (we assume that each genotype class of subjects have the same standard 
                            deviation). This level of noise can be based on estimates from previous study or pilot data. "),

                    tags$br(),
                    tags$br(),

                    paste0("More details can be found "),
                    HTML("<a target='_blank' onclick=","customHref('about')>here</a>"),
                    paste0("and in the online manual of the function ‘powerEQTL.ANOVA’ in R package "),
                    tags$a( "‘powerEQTL’", href = "https://CRAN.R-project.org/package=powerEQTL", target="_blank"),
                    
                    tags$br(),
                )
            })
        }
        else
        {
            updateNumericInput(session, inputId = "delta1", label = "Slope of the regression line (β1)", value = input$delta1)
            updateNumericInput(session, inputId = "sigma1", 
                               label = "Standard deviation of gene expression (σy)",
                               value = input$sigma1, min = 0.001, max = 1, step = 0.001)
            hide("delta2")
            output$description2 <- renderUI({
                tags$div(
                    
                    withMathJax(),
                    
                    # tags$br(),
                    
                    paste0("The figure above shows the statistical power of an eQTL analysis as a function of the minor
                            allele frequency (MAF) with different sample sizes.  Simple linear regression model was used 
                            to test the linear association between gene expression and genotype.
                            In this example, we model the expression data as log-normally distributed with a log standard
                            deviation of ", input$sigma1, ". This level of noise can be based on estimates from previous
                            study or pilot data. "),

                    tags$br(),
                    tags$br(),

                    paste0("More details can be found "),
                    HTML("<a onclick=","customHref('about') formtarget='_blank'>here</a>"),
                    paste0("and in the online manual of the function ‘powerEQTL.SLR’ in R package "),
                    tags$a( "‘powerEQTL’", href = "https://CRAN.R-project.org/package=powerEQTL", target="_blank"),
                    tags$br(),
                    )
            })
        }
        ## update plot
        output$tissue <- renderPlot(tissue())
    })
    ## Interactive plot for tissue
    observeEvent(input$t_plot_hover$y,{
        # get hover coordinate
        power_t <<- input$t_plot_hover$y
        
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
            
        }
        
        ## update plot
        output$tissue <- renderPlot({
                # set up graph
                xrange <- range(mafVec1()*100)
                yrange <- c(0:1)
                colors <- hcl.colors(sizeLen.data(),"Set 2")
                par(mar = c(5, 5, 1, 2), cex.lab = 1.5)
                
                ## 	 # treating each cell independantly 
                plot(xrange, yrange, xlim=range(xrange), log='x', type="n",
                     xlab="MAF (%)",
                     ylab="Power",
                     xaxt = "n", yaxt = "n")
                
                axis(1, cex.axis = 1.5)
                axis(2, cex.axis = 1.5)
                
                abline(v=0, h=seq(0,1,.1), lty=2, col="grey89")
                abline(h=0, v=c(1:10), lty=2,col="grey89")
                abline(h=power_t, lty=2, lwd = 1.5, col = 1)
                text(mafVec1()[1]*100, power_t,labels = as.character(round(power_t,2)), cex=1.5, adj=0.5, pos = 1, offset = -1)  #power level
                
                # add power curves
                for (i in 1:sizeLen.data()){
                    lines(mafVec1()*100, power_tissue[i,], type="l", lwd=4, col=colors[i])
                    pos1=which(power_tissue[i,]>=power_t)[1]
                    abline(v=mafVec1()[pos1]*100, col=colors[i])
                    text(mafVec1()[pos1]*100,power_t,labels = as.character(round(mafVec1()[pos1]*100,2)), srt = 90,cex=1.5, adj=0.5)
                }
                
                legend(input$pos_t, title="Sample size (n)", legend = sizeVec.data(),
                       title.col='black',text.col=colors,cex =1.5, bty='n', pch = 15, col = colors)
            })

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
            
            mtext(paste("linear mixed effects model (FWER=",input$alpha, ", nTests=", input$nTest,", slope=",input$delta,", m="
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
            
            legend(input$pos, title="Sample size (n)", legend = sizeVec1.data(),
                   title.col='black',text.col=colors,cex =1, bty='n', pch = 15, col = colors)
            
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
                mtext(paste("One-way unbalanced ANOVA (FWER=",input$alpha1, ", nTests=", input$nTest1,", delta=",input$delta1,", sigma=", input$sigma1 ,")", sep = ""),3)
                
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
                mtext(paste("Simple Linear Regression (FWER=",input$alpha1, ", nTests=", input$nTest1,", beta=",input$delta1,", sigma=", input$sigma1 ,")", sep = ""),3)
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
            legend(input$pos_t, title="Sample size (n)", legend = sizeVec.data(),
                   title.col='black',text.col=colors,cex =1, bty='n', pch = 15, col = colors)
            
            dev.off()
        }
    )
    
    
    ## Single cell Plot Explanations
    output$description <- renderUI({
        tags$div(
        
        tags$br(),
        
        HTML(paste0("The figure above shows the statistical power of an eQTL analysis as a function of the minor allele
                    frequency (MAF) with different sample sizes. Simple linear mixed effects regression model was used 
                    to characterize the linear association between gene expression and genotype. We assume expression 
                    levels are normally distributed after appropriate pre-processing. In this example, We model the 
                    expression data as log-normally distributed with a log standard deviation of ", input$sigma, " 
                    and an intra-class correlation of ", input$rho, ". This level of noise can be based on estimates 
                    from previous study or pilot data. ")),
                    
        tags$br(),
        tags$br(),  
        
        paste0("More details can be found "),
        HTML("<a target='_blank' onclick=","customHref('about')>here</a>"),
        paste0("and in the online manual of the function ‘powerEQTL.scRNAseq’ in R package "),
        tags$a( "‘powerEQTL’", href = "https://CRAN.R-project.org/package=powerEQTL", target="_blank"),
        tags$br()
        )
    })
    
    ## Power calculator for single-cell eQTL
    ## hidden input options
    observeEvent(input$btn_scest, {
        
        if (input$btn_scest %% 2 == 1) {
            txt <- "Less Options"
            show("rho_est")
            show("rho_est_t")
            show("m_est")
            show("m_est_t")
            show("sigma_est")
            show("sigma_est_t")
            show("FWER_est")
            show("FWER_est_t")
            show("nTest_est")
            show("nTest_est_t")
        } 
        else {
            txt <- "More Options"
            hide("rho_est")
            hide("rho_est_t")
            hide("m_est")
            hide("m_est_t")
            hide("sigma_est")
            hide("sigma_est_t")
            hide("FWER_est")
            hide("FWER_est_t")
            hide("nTest_est")
            hide("nTest_est_t")
        }
        updateActionButton(session, "btn_scest", txt)
    })
    ## Approximate sample sizes with given inputs
    observeEvent(input$power_est,{
        max = powerEQTL.scRNAseq(
            slope = input$slope_est, 
            n = NULL,
            m = input$m_est, 
            power = 0.99999,
            sigma.y = input$sigma_est, 
            MAF = input$maf_est, 
            rho = input$rho_est, 
            FWER = input$FWER_est,
            nTests = input$nTest_est)
        n = powerEQTL.scRNAseq(
            slope = input$slope_est, 
            n = NULL,
            m = input$m_est, 
            power = input$power_est,
            sigma.y = input$sigma_est, 
            MAF = input$maf_est, 
            rho = input$rho_est, 
            FWER = input$FWER_est,
            nTests = input$nTest_est)
        updateSliderInput(session, inputId = "n_est", label = "Number of subjects needed",
                          value = as.numeric(n), max = ceiling(as.numeric(max)))
    })
    observeEvent(input$n_est,{
        power = powerEQTL.scRNAseq(
            power = NULL,
            slope = input$slope_est,
            m = input$m_est,
            n = input$n_est,
            sigma.y = input$sigma_est,
            MAF = input$maf_est,
            rho = input$rho_est,
            FWER = input$FWER_est,
            nTests = input$nTest_est)
        updateSliderInput(session, inputId = "power_est", label = "Desired power level",
                          value = as.numeric(power), min = 0.001, max = 1)
    })
    observeEvent(input$maf_est,{
        max = powerEQTL.scRNAseq(
            slope = input$slope_est, 
            n = NULL,
            m = input$m_est, 
            power = 0.99999,
            sigma.y = input$sigma_est, 
            MAF = input$maf_est, 
            rho = input$rho_est, 
            FWER = input$FWER_est,
            nTests = input$nTest_est)
        n = powerEQTL.scRNAseq(
            slope = input$slope_est,
            n = NULL,
            m = input$m_est,
            power = input$power_est,
            sigma.y = input$sigma_est,
            MAF = input$maf_est,
            rho = input$rho_est,
            FWER = input$FWER_est,
            nTests = input$nTest_est)
        updateSliderInput(session, inputId = "n_est", label = "Number of subjects needed",
                          value = as.numeric(n), max = ceiling(as.numeric(max)))
    })
    observeEvent(input$slope_est,{
        max = powerEQTL.scRNAseq(
            slope = input$slope_est, 
            n = NULL,
            m = input$m_est, 
            power = 0.99999,
            sigma.y = input$sigma_est, 
            MAF = input$maf_est, 
            rho = input$rho_est, 
            FWER = input$FWER_est,
            nTests = input$nTest_est)
        n = powerEQTL.scRNAseq(
            slope = input$slope_est,
            n = NULL,
            m = input$m_est,
            power = input$power_est,
            sigma.y = input$sigma_est,
            MAF = input$maf_est,
            rho = input$rho_est,
            FWER = input$FWER_est,
            nTests = input$nTest_est)
        updateSliderInput(session, inputId = "n_est", label = "Number of subjects needed",
                          value = as.numeric(n), max = ceiling(as.numeric(max)))
        updateSliderInput(session, inputId = "slope_est", label = "Slope of the regression line (β1)",
                          value = input$slope_est, min = 0.01, max = max(1, input$slope_est*2))
    })
    observeEvent(input$sigma_est,{
        max = powerEQTL.scRNAseq(
            slope = input$slope_est, 
            n = NULL,
            m = input$m_est, 
            power = 0.99999,
            sigma.y = input$sigma_est, 
            MAF = input$maf_est, 
            rho = input$rho_est, 
            FWER = input$FWER_est,
            nTests = input$nTest_est)
        n = powerEQTL.scRNAseq(
            slope = input$slope_est,
            n = NULL,
            m = input$m_est,
            power = input$power_est,
            sigma.y = input$sigma_est,
            MAF = input$maf_est,
            rho = input$rho_est,
            FWER = input$FWER_est,
            nTests = input$nTest_est)
        updateSliderInput(session, inputId = "n_est", label = "Number of subjects needed",
                          value = as.numeric(n), max = ceiling(as.numeric(max)))
        updateSliderInput(session, inputId = "sigma_est", 
                          value = input$sigma_est, min = 0.01, max = max(10, input$sigma_est*2))
    })
    observeEvent(input$FWER_est,{
        max = powerEQTL.scRNAseq(
            slope = input$slope_est, 
            n = NULL,
            m = input$m_est, 
            power = 0.99999,
            sigma.y = input$sigma_est, 
            MAF = input$maf_est, 
            rho = input$rho_est, 
            FWER = input$FWER_est,
            nTests = input$nTest_est)
        n = powerEQTL.scRNAseq(
            slope = input$slope_est,
            n = NULL,
            m = input$m_est,
            power = input$power_est,
            sigma.y = input$sigma_est,
            MAF = input$maf_est,
            rho = input$rho_est,
            FWER = input$FWER_est,
            nTests = input$nTest_est)
        updateSliderInput(session, inputId = "n_est", label = "Number of subjects needed",
                          value = as.numeric(n), max = ceiling(as.numeric(max)))
    })
    observeEvent(input$nTest_est,{
        max = powerEQTL.scRNAseq(
            slope = input$slope_est, 
            n = NULL,
            m = input$m_est, 
            power = 0.99999,
            sigma.y = input$sigma_est, 
            MAF = input$maf_est, 
            rho = input$rho_est, 
            FWER = input$FWER_est,
            nTests = input$nTest_est)
        n = powerEQTL.scRNAseq(
            slope = input$slope_est,
            n = NULL,
            m = input$m_est,
            power = input$power_est,
            sigma.y = input$sigma_est,
            MAF = input$maf_est,
            rho = input$rho_est,
            FWER = input$FWER_est,
            nTests = input$nTest_est)
        updateSliderInput(session, inputId = "n_est", label = "Number of subjects needed",
                          value = as.numeric(n), max = ceiling(as.numeric(max)))
        updateSliderInput(session, inputId = "nTest_est", label = "Total number of tests (nTests)",
                          value = input$nTest_est, min = 0, max = max(10e7, input$nTest_est*2))
    })
    observeEvent(input$rho_est,{
        max = powerEQTL.scRNAseq(
            slope = input$slope_est, 
            n = NULL,
            m = input$m_est, 
            power = 0.99999,
            sigma.y = input$sigma_est, 
            MAF = input$maf_est, 
            rho = input$rho_est, 
            FWER = input$FWER_est,
            nTests = input$nTest_est)
        n = powerEQTL.scRNAseq(
            slope = input$slope_est,
            n = NULL,
            m = input$m_est,
            power = input$power_est,
            sigma.y = input$sigma_est,
            MAF = input$maf_est,
            rho = input$rho_est,
            FWER = input$FWER_est,
            nTests = input$nTest_est)
        updateSliderInput(session, inputId = "n_est", label = "Number of subjects needed",
                          value = as.numeric(n), max = ceiling(as.numeric(max)))
    })
    observeEvent(input$m_est,{
        max = powerEQTL.scRNAseq(
            slope = input$slope_est, 
            n = NULL,
            m = input$m_est, 
            power = 0.99999,
            sigma.y = input$sigma_est, 
            MAF = input$maf_est, 
            rho = input$rho_est, 
            FWER = input$FWER_est,
            nTests = input$nTest_est)
        n = powerEQTL.scRNAseq(
            slope = input$slope_est,
            n = NULL,
            m = input$m_est,
            power = input$power_est,
            sigma.y = input$sigma_est,
            MAF = input$maf_est,
            rho = input$rho_est,
            FWER = input$FWER_est,
            nTests = input$nTest_est)
        updateSliderInput(session, inputId = "n_est", label = "Number of subjects needed",
                          value = as.numeric(n), max = ceiling(as.numeric(max)))
        updateSliderInput(session, inputId = "m_est", 
                          value = input$m_est, min = 500, max = max(100000, input$m_est*2))
    })
    
    ## Power calculator for tissue eQTL

    # Reactive Input Title for delta for different models
    observeEvent(input$radio_test,{
        if(input$radio_test == "One-way unbalanced ANOVA")
        {
            hide("slope_test")
            hide("slope_test_t")
            if (input$btn_test %% 2 == 1)
            {
                show("delta1_test")
                show("delta1_test_t")
                show("delta2_test")
                show("delta2_test_t")
            }
            max = powerEQTL.ANOVA(
                n = NULL,
                deltaVec = c(input$delta1_test, input$delta2_test),
                power = 0.999,
                sigma = input$sigma_test,
                MAF = input$maf_test,
                FWER = input$FWER_test,
                nTests = input$nTest_test)
            n = powerEQTL.ANOVA(
                n = NULL,
                deltaVec = c(input$delta1_test, input$delta2_test),
                power = input$power_test,
                sigma = input$sigma_test,
                MAF = input$maf_test,
                FWER = input$FWER_test,
                nTests = input$nTest_test)
            
            updateSliderInput(session, inputId = "n_test", label = "Num subjects needed",
                              value = as.numeric(n), max = ceiling(as.numeric(max)))
            updateSliderInput(session, inputId = "sigma_test", 
                               label = "Standard deviation of random error (σ)",
                               value = input$sigma_test, min = 0.001, max = 1, step = 0.001)
            output$Explanation1 <- renderUI({
                tags$div(
                    tags$br(),
                    paste0("* - Minimum number of subjects needed for designated power level."),
                    tags$br(),
                    tags$br(),
                    paste0("More details can be found "),
                    HTML("<a onclick=","customHref('about'), target='_blank'>here</a>"),
                    paste0("and in the online manual of the function 'powerEQTL.ANOVA' in R package "),
                    tags$a( "‘powerEQTL’", href = "https://CRAN.R-project.org/package=powerEQTL", target="_blank"),
                    tags$br(),
                    tags$br(),
                )
            })
        }
        else
        {
            show("slope_test")
            show("slope_test_t")
            hide("delta1_test")
            hide("delta1_test_t")
            hide("delta2_test")
            hide("delta2_test_t")
            max = powerEQTL.SLR(
                n = NULL,
                slope = input$slope_test,
                power = 0.999,
                sigma.y = input$sigma_test,
                MAF = input$maf_test,
                FWER = input$FWER_test,
                nTests = input$nTest_test)
            n = powerEQTL.SLR(
                n = NULL,
                slope = input$slope_test,
                power = input$power_test,
                sigma.y = input$sigma_test,
                MAF = input$maf_test,
                FWER = input$FWER_test,
                nTests = input$nTest_test)
            updateSliderInput(session, inputId = "n_test", label = "Num subjects needed",
                              value = as.numeric(n), max = ceiling(as.numeric(max)))
            updateNumericInput(session, inputId = "sigma_test", 
                               label = "Standard deviation of gene expression (σy)",
                               value = input$sigma_test, min = 0.001, max = 1, step = 0.001)
            
            output$Explanation1 <- renderUI({
                tags$div(
                    tags$br(),
                    paste0("* - Minimum number of subjects needed for designated power level."),
                    tags$br(),
                    tags$br(),
                    paste0("More details can be found "),
                    HTML("<a onclick=","customHref('about'), target='_blank'>here</a>"),
                    paste0("and in the online manual of the function 'powerEQTL.SLR' in R package "),
                    tags$a( "‘powerEQTL’", href = "https://CRAN.R-project.org/package=powerEQTL", target="_blank"),
                    tags$br(),
                    tags$br(),
                )
            })
            
            
        }

    })
    
    ## hidden input options
    observeEvent(input$btn_test, {
        
        if (input$btn_test %% 2 == 1) {
            txt <- "Less Options"
            if(input$radio_test == "One-way unbalanced ANOVA")
            {
                show("delta1_test")
                show("delta1_test_t")
                show("delta2_test")
                show("delta2_test_t")
            }
            
            show("sigma_test")
            show("sigma_test_t")
            show("FWER_test")
            show("FWER_test_t")
            show("nTest_test")
            show("nTest_test_t")
            
        } 
        else {
            txt <- "More Options"
            if(input$radio_test == "One-way unbalanced ANOVA")
            {
                toggle("delta1_test")
                toggle("delta1_test_t")
                toggle("delta2_test")
                toggle("delta2_test_t")
            }
            
            hide("sigma_test")
            hide("sigma_test_t")
            hide("FWER_test")
            hide("FWER_test_t")
            hide("nTest_test")
            hide("nTest_test_t")

        }
        updateActionButton(session, "btn_test", txt)
    })

    ## Approximate sample sizes with given inputs
    observeEvent(input$power_test,{
        if(input$radio_test == "One-way unbalanced ANOVA")
        {
            max = powerEQTL.ANOVA(
                n = NULL,
                deltaVec = c(input$delta1_test, input$delta2_test),
                power = 0.999,
                sigma = input$sigma_test,
                MAF = input$maf_test,
                FWER = input$FWER_test,
                nTests = input$nTest_test)
            n = powerEQTL.ANOVA(
                n = NULL,
                deltaVec = c(input$delta1_test, input$delta2_test),
                power = input$power_test,
                sigma = input$sigma_test,
                MAF = input$maf_test,
                FWER = input$FWER_test,
                nTests = input$nTest_test)
            
            updateSliderInput(session, inputId = "n_test", label = "Num subjects needed",
                              value = as.numeric(n), max = ceiling(as.numeric(max)))
        }
        else
        {
            max = powerEQTL.SLR(
                n = NULL,
                slope = input$slope_test,
                power = 0.999,
                sigma.y = input$sigma_test,
                MAF = input$maf_test,
                FWER = input$FWER_test,
                nTests = input$nTest_test)
            n = powerEQTL.SLR(
                n = NULL,
                slope = input$slope_test,
                power = input$power_test,
                sigma.y = input$sigma_test,
                MAF = input$maf_test,
                FWER = input$FWER_test,
                nTests = input$nTest_test)
            updateSliderInput(session, inputId = "n_test", label = "Num subjects needed",
                              value = as.numeric(n), max = ceiling(as.numeric(max)))
        }
    })
    observeEvent(input$maf_test,{
        if(input$radio_test == "One-way unbalanced ANOVA")
        {
            max = powerEQTL.ANOVA(
                n = NULL,
                deltaVec = c(input$delta1_test, input$delta2_test),
                power = 0.999,
                sigma = input$sigma_test,
                MAF = input$maf_test,
                FWER = input$FWER_test,
                nTests = input$nTest_test)
            n = powerEQTL.ANOVA(
                n = NULL,
                deltaVec = c(input$delta1_test, input$delta2_test),
                power = input$power_test,
                sigma = input$sigma_test,
                MAF = input$maf_test,
                FWER = input$FWER_test,
                nTests = input$nTest_test)
            
            updateSliderInput(session, inputId = "n_test", label = "Num subjects needed",
                              value = ceiling(as.numeric(n)), max = ceiling(as.numeric(max)))
        }
        else
        {
            max = powerEQTL.SLR(
                n = NULL,
                slope = input$slope_test,
                power = 0.999,
                sigma.y = input$sigma_test,
                MAF = input$maf_test,
                FWER = input$FWER_test,
                nTests = input$nTest_test)
            n = powerEQTL.SLR(
                n = NULL,
                slope = input$slope_test,
                power = input$power_test,
                sigma.y = input$sigma_test, 
                MAF = input$maf_test, 
                FWER = input$FWER_test,
                nTests = input$nTest_test)
            updateSliderInput(session, inputId = "n_test", label = "Num subjects needed",
                              value = ceiling(as.numeric(n)), max = ceiling(as.numeric(max)))
        }
    })
    observeEvent(input$slope_test,{
        if(input$radio_test != "One-way unbalanced ANOVA")
        {
            max = powerEQTL.SLR(
                n = NULL,
                slope = input$slope_test,
                power = 0.999,
                sigma.y = input$sigma_test,
                MAF = input$maf_test,
                FWER = input$FWER_test,
                nTests = input$nTest_test)
            n = powerEQTL.SLR(
                n = NULL,
                slope = input$slope_test,
                power = input$power_test,
                sigma.y = input$sigma_test, 
                MAF = input$maf_test, 
                FWER = input$FWER_test,
                nTests = input$nTest_test)
            updateSliderInput(session, inputId = "n_test", label = "Num subjects needed",
                              value = ceiling(as.numeric(n)), max = ceiling(as.numeric(max)))
            updateSliderInput(session, inputId = "slope_test", label = "Slope of the regression line (β1)",
                        value = input$slope_test, min = 0.01, max = max(1, input$slope_test*2))
        }
    })
    observeEvent(input$n_test,{
        if(input$radio_test == "One-way unbalanced ANOVA")
        {
            p = powerEQTL.ANOVA(
                n = input$n_test,
                deltaVec = c(input$delta1_test, input$delta2_test),
                power = NULL,
                sigma = input$sigma_test,
                MAF = input$maf_test,
                FWER = input$FWER_test,
                nTests = input$nTest_test)
            updateSliderInput(session, inputId = "power_test", label = "Desired power level",
                              value = round(as.numeric(p),2))
        }
        else
        {
            p = powerEQTL.SLR(
                n = input$n_test,
                slope = input$slope_test,
                power = NULL,
                sigma.y = input$sigma_test, 
                MAF = input$maf_test, 
                FWER = input$FWER_test,
                nTests = input$nTest_test)
            updateSliderInput(session, inputId = "power_test", label = "Desired power level",
                              value = round(as.numeric(p),2))
        }
        
    })
    observeEvent(input$delta1_test,{
        if(input$radio_test == "One-way unbalanced ANOVA")
        {
            max = powerEQTL.ANOVA(
                n = NULL,
                deltaVec = c(input$delta1_test, input$delta2_test),
                power = 0.999,
                sigma = input$sigma_test,
                MAF = input$maf_test,
                FWER = input$FWER_test,
                nTests = input$nTest_test)
            n = powerEQTL.ANOVA(
                n = NULL,
                deltaVec = c(input$delta1_test, input$delta2_test),
                power = input$power_test,
                sigma = input$sigma_test,
                MAF = input$maf_test,
                FWER = input$FWER_test,
                nTests = input$nTest_test)
            
            updateSliderInput(session, inputId = "n_test", label = "Num subjects needed",
                              value = ceiling(as.numeric(n)), max = ceiling(as.numeric(max)))
            updateSliderInput(session, inputId = "delta1_test", label = "Mean difference of gene expression (δ1)",
                        value = input$delta1_test, min = 0.01, max = max(10, input$delta1_test*2))
        }
    })
    observeEvent(input$delta2_test,{
        if(input$radio_test == "One-way unbalanced ANOVA")
        {
            max = powerEQTL.ANOVA(
                n = NULL,
                deltaVec = c(input$delta1_test, input$delta2_test),
                power = 0.999,
                sigma = input$sigma_test,
                MAF = input$maf_test,
                FWER = input$FWER_test,
                nTests = input$nTest_test)
            n = powerEQTL.ANOVA(
                n = NULL,
                deltaVec = c(input$delta1_test, input$delta2_test),
                power = input$power_test,
                sigma = input$sigma_test,
                MAF = input$maf_test,
                FWER = input$FWER_test,
                nTests = input$nTest_test)
            
            updateSliderInput(session, inputId = "n_test", label = "Num subjects needed",
                              value = ceiling(as.numeric(n)), max = ceiling(as.numeric(max)))
            updateSliderInput(session, inputId = "delta2_test", label = "Mean difference of gene expression (δ2)",
                              value = input$delta2_test, min = 0.01, max = max(10, input$delta2_test*2))
        }
    })
    observeEvent(input$sigma_test,{
        if(input$radio_test == "One-way unbalanced ANOVA")
        {
            max = powerEQTL.ANOVA(
                n = NULL,
                deltaVec = c(input$delta1_test, input$delta2_test),
                power = 0.999,
                sigma = input$sigma_test,
                MAF = input$maf_test,
                FWER = input$FWER_test,
                nTests = input$nTest_test)
            n = powerEQTL.ANOVA(
                n = NULL,
                deltaVec = c(input$delta1_test, input$delta2_test),
                power = input$power_test,
                sigma = input$sigma_test,
                MAF = input$maf_test,
                FWER = input$FWER_test,
                nTests = input$nTest_test)
            
            updateSliderInput(session, inputId = "n_test", label = "Num ubjects needed",
                              value = ceiling(as.numeric(n)), max = ceiling(as.numeric(max)))
            updateSliderInput(session, inputId = "sigma_test", 
                              value = input$sigma_test, min = 0.01, max = max(10, input$sigma_test*2))
        }
        else
        {
            max = powerEQTL.SLR(
                n = NULL,
                slope = input$slope_test,
                power = 0.999,
                sigma.y = input$sigma_test,
                MAF = input$maf_test,
                FWER = input$FWER_test,
                nTests = input$nTest_test)
            n = powerEQTL.SLR(
                n = NULL,
                slope = input$slope_test,
                power = input$power_test,
                sigma.y = input$sigma_test, 
                MAF = input$maf_test, 
                FWER = input$FWER_test,
                nTests = input$nTest_test)
            updateSliderInput(session, inputId = "n_test", label = "Num subjects needed",
                              value = ceiling(as.numeric(n)), max = ceiling(as.numeric(max)))
            updateSliderInput(session, inputId = "sigma_test", 
                              value = input$sigma_test, min = 0.01, max = max(10, input$sigma_test*2))
        }
    })
    observeEvent(input$FWER_test,{
        if(input$radio_test == "One-way unbalanced ANOVA")
        {
            max = powerEQTL.ANOVA(
                n = NULL,
                deltaVec = c(input$delta1_test, input$delta2_test),
                power = 0.999,
                sigma = input$sigma_test,
                MAF = input$maf_test,
                FWER = input$FWER_test,
                nTests = input$nTest_test)
            n = powerEQTL.ANOVA(
                n = NULL,
                deltaVec = c(input$delta1_test, input$delta2_test),
                power = input$power_test,
                sigma = input$sigma_test,
                MAF = input$maf_test,
                FWER = input$FWER_test,
                nTests = input$nTest_test)
            
            updateSliderInput(session, inputId = "n_test", label = "Num subjects needed",
                              value = ceiling(as.numeric(n)), max = ceiling(as.numeric(max)))
        }
        else
        {
            max = powerEQTL.SLR(
                n = NULL,
                slope = input$slope_test,
                power = 0.999,
                sigma.y = input$sigma_test,
                MAF = input$maf_test,
                FWER = input$FWER_test,
                nTests = input$nTest_test)
            n = powerEQTL.SLR(
                n = NULL,
                slope = input$slope_test,
                power = input$power_test,
                sigma.y = input$sigma_test, 
                MAF = input$maf_test, 
                FWER = input$FWER_test,
                nTests = input$nTest_test)
            updateSliderInput(session, inputId = "n_test", label = "Num subjects needed",
                              value = ceiling(as.numeric(n)), max = ceiling(as.numeric(max)))
        }
    })
    observeEvent(input$nTest_test,{
        if(input$radio_test == "One-way unbalanced ANOVA")
        {
            max = powerEQTL.ANOVA(
                n = NULL,
                deltaVec = c(input$delta1_test, input$delta2_test),
                power = 0.999,
                sigma = input$sigma_test,
                MAF = input$maf_test,
                FWER = input$FWER_test,
                nTests = input$nTest_test)
            n = powerEQTL.ANOVA(
                n = NULL,
                deltaVec = c(input$delta1_test, input$delta2_test),
                power = input$power_test,
                sigma = input$sigma_test,
                MAF = input$maf_test,
                FWER = input$FWER_test,
                nTests = input$nTest_test)
            
            updateSliderInput(session, inputId = "n_test", label = "Num subjects needed",
                              value = ceiling(as.numeric(n)), max = ceiling(as.numeric(max)))
            updateSliderInput(session, inputId = "nTest_test", label = "Total number of tests (nTests)",
                              value = input$nTest_test, min = 0, max = max(10e7, input$nTest_test*2))
        }
        else
        {
            max = powerEQTL.SLR(
                n = NULL,
                slope = input$slope_test,
                power = 0.999,
                sigma.y = input$sigma_test,
                MAF = input$maf_test,
                FWER = input$FWER_test,
                nTests = input$nTest_test)
            n = powerEQTL.SLR(
                n = NULL,
                slope = input$slope_test,
                power = input$power_test,
                sigma.y = input$sigma_test, 
                MAF = input$maf_test, 
                FWER = input$FWER_test,
                nTests = input$nTest_test)
            updateSliderInput(session, inputId = "n_test", label = "Num subjects needed",
                              value = ceiling(as.numeric(n)), max = ceiling(as.numeric(max)))
            updateSliderInput(session, inputId = "nTest_test", label = "Total number of tests (nTests)",
                              value = input$nTest_test, min = 0, max = max(10e7, input$nTest_test*2))
        }
    })
    
    ## Update slider and text inputs accordingly
    # Tissue
    observeEvent(input$power_test_t,{
        if(as.numeric(input$power_test_t) != input$power_test)
        {
            updateSliderInput(
                session = session,
                inputId = 'power_test',
                value = input$power_test_t
            )
        }
        
        
    })
    observeEvent(input$power_test,{
        if(as.numeric(input$power_test_t) != input$power_test)
        {
            updateTextInput(
                session = session,
                inputId = 'power_test_t',
                value = input$power_test
            )
            
        }
        
    })
    observeEvent(input$n_test_t,{
        if(as.numeric(input$n_test_t) != input$n_test)
        {
            updateSliderInput(
                session = session,
                inputId = 'n_test',
                value = input$n_test_t
            )
        }
        
        
    })
    observeEvent(input$n_test,{
        if(as.numeric(input$n_test_t) != input$n_test)
        {
            updateTextInput(
                session = session,
                inputId = 'n_test_t',
                value = input$n_test
            )
            
        }
        
    })
    observeEvent(input$maf_test_t,{
        if(as.numeric(input$maf_test_t) != input$maf_test)
        {
            updateSliderInput(
                session = session,
                inputId = 'maf_test',
                value = input$maf_test_t
            )
        }
    })
    observeEvent(input$maf_test,{
        if(as.numeric(input$maf_test_t) != input$maf_test)
        {
            updateTextInput(
                session = session,
                inputId = 'maf_test_t',
                value = input$maf_test
            )
            
        }
        
    })
    observeEvent(input$slope_test_t,{
        if(as.numeric(input$slope_test_t) != input$slope_test)
        {
            updateSliderInput(
                session = session,
                inputId = 'slope_test',
                value = input$slope_test_t
            )
        }
    })
    observeEvent(input$slope_test,{
        if(as.numeric(input$slope_test_t) != input$slope_test)
        {
            updateTextInput(
                session = session,
                inputId = 'slope_test_t',
                value = input$slope_test
            )
            
        }
        
    })
    observeEvent(input$delta1_test_t,{
        if(as.numeric(input$delta1_test_t) != input$delta1_test)
        {
            updateSliderInput(
                session = session,
                inputId = 'delta1_test',
                value = input$delta1_test_t
            )
        }
    })
    observeEvent(input$delta1_test,{
        if(as.numeric(input$delta1_test_t) != input$delta1_test)
        {
            updateTextInput(
                session = session,
                inputId = 'delta1_test_t',
                value = input$delta1_test
            )
            
        }
        
    })
    observeEvent(input$delta2_test_t,{
        if(as.numeric(input$delta2_test_t) != input$delta2_test)
        {
            updateSliderInput(
                session = session,
                inputId = 'delta2_test',
                value = input$delta2_test_t
            )
        }
    })
    observeEvent(input$delta2_test,{
        if(as.numeric(input$delta2_test_t) != input$delta2_test)
        {
            updateTextInput(
                session = session,
                inputId = 'delta2_test_t',
                value = input$delta2_test
            )
            
        }
    })
    observeEvent(input$sigma_test_t,{
        if(as.numeric(input$sigma_test_t) != input$sigma_test)
        {
            updateSliderInput(
                session = session,
                inputId = 'sigma_test',
                value = input$sigma_test_t
            )
        }
    })
    observeEvent(input$sigma_test,{
        if(as.numeric(input$sigma_test_t) != input$sigma_test)
        {
            updateTextInput(
                session = session,
                inputId = 'sigma_test_t',
                value = input$sigma_test
            )
            
        }
    })
    observeEvent(input$FWER_test_t,{
        if(as.numeric(input$FWER_test_t) != input$FWER_test)
        {
            updateSliderInput(
                session = session,
                inputId = 'FWER_test',
                value = input$FWER_test_t
            )
        }
    })
    observeEvent(input$FWER_test,{
        if(as.numeric(input$FWER_test_t) != input$FWER_test)
        {
            updateTextInput(
                session = session,
                inputId = 'FWER_test_t',
                value = input$FWER_test
            )
            
        }
    })
    observeEvent(input$nTest_test_t,{
        if(as.numeric(input$nTest_test_t) != input$nTest_test)
        {
            updateSliderInput(
                session = session,
                inputId = 'nTest_test',
                value = input$nTest_test_t
            )
        }
    })
    observeEvent(input$nTest_test,{
        if(as.numeric(input$nTest_test_t) != input$nTest_test)
        {
            updateTextInput(
                session = session,
                inputId = 'nTest_test_t',
                value = input$nTest_test
            )
            
        }
    })
    # Single cell
    observeEvent(input$power_est_t,{
        if(as.numeric(input$power_est_t) != input$power_est)
        {
            updateSliderInput(
                session = session,
                inputId = 'power_est',
                value = input$power_est_t
            )
        }
        
        
    })
    observeEvent(input$power_est,{
        if(as.numeric(input$power_est_t) != input$power_est)
        {
            updateTextInput(
                session = session,
                inputId = 'power_est_t',
                value = input$power_est
            )
            
        }
        
    })
    observeEvent(input$n_est_t,{
        if(as.numeric(input$n_est_t) != input$n_est)
        {
            updateSliderInput(
                session = session,
                inputId = 'n_est',
                value = input$n_est_t
            )
        }
        
        
    })
    observeEvent(input$n_est,{
        if(as.numeric(input$n_est_t) != input$n_est)
        {
            updateTextInput(
                session = session,
                inputId = 'n_est_t',
                value = input$n_est
            )
            
        }
        
    })
    observeEvent(input$maf_est_t,{
        if(as.numeric(input$maf_est_t) != input$maf_est)
        {
            updateSliderInput(
                session = session,
                inputId = 'maf_est',
                value = input$maf_est_t
            )
        }
    })
    observeEvent(input$maf_est,{
        if(as.numeric(input$maf_est_t) != input$maf_est)
        {
            updateTextInput(
                session = session,
                inputId = 'maf_est_t',
                value = input$maf_est
            )
            
        }
        
    })
    observeEvent(input$slope_est_t,{
        if(as.numeric(input$slope_est_t) != input$slope_est)
        {
            updateSliderInput(
                session = session,
                inputId = 'slope_est',
                value = input$slope_est_t
            )
        }
    })
    observeEvent(input$slope_est,{
        if(as.numeric(input$slope_est_t) != input$slope_est)
        {
            updateTextInput(
                session = session,
                inputId = 'slope_est_t',
                value = input$slope_est
            )
            
        }
        
    })
    observeEvent(input$sigma_est_t,{
        if(as.numeric(input$sigma_est_t) != input$sigma_est)
        {
            updateSliderInput(
                session = session,
                inputId = 'sigma_est',
                value = input$sigma_est_t
            )
        }
    })
    observeEvent(input$sigma_est,{
        if(as.numeric(input$sigma_est_t) != input$sigma_est)
        {
            updateTextInput(
                session = session,
                inputId = 'sigma_est_t',
                value = input$sigma_est
            )
            
        }
    })
    observeEvent(input$FWER_est_t,{
        if(as.numeric(input$FWER_est_t) != input$FWER_est)
        {
            updateSliderInput(
                session = session,
                inputId = 'FWER_est',
                value = input$FWER_est_t
            )
        }
    })
    observeEvent(input$FWER_est,{
        if(as.numeric(input$FWER_est_t) != input$FWER_est)
        {
            updateTextInput(
                session = session,
                inputId = 'FWER_est_t',
                value = input$FWER_est
            )
            
        }
    })
    observeEvent(input$nTest_est_t,{
        if(as.numeric(input$nTest_est_t) != input$nTest_est)
        {
            updateSliderInput(
                session = session,
                inputId = 'nTest_est',
                value = input$nTest_est_t
            )
        }
    })
    observeEvent(input$nTest_est,{
        if(as.numeric(input$nTest_est_t) != input$nTest_est)
        {
            updateTextInput(
                session = session,
                inputId = 'nTest_est_t',
                value = input$nTest_est
            )
            
        }
    })
    observeEvent(input$rho_est_t,{
        if(as.numeric(input$rho_est_t) != input$rho_est)
        {
            updateSliderInput(
                session = session,
                inputId = 'rho_est',
                value = input$rho_est_t
            )
        }
    })
    observeEvent(input$rho_est,{
        if(as.numeric(input$rho_est_t) != input$rho_est)
        {
            updateTextInput(
                session = session,
                inputId = 'rho_est_t',
                value = input$rho_est
            )
            
        }
    })
    observeEvent(input$m_est_t,{
        if(as.numeric(input$m_est_t) != input$m_est)
        {
            updateSliderInput(
                session = session,
                inputId = 'm_est',
                value = input$m_est_t
            )
        }
    })
    observeEvent(input$m_est,{
        if(as.numeric(input$m_est_t) != input$m_est)
        {
            updateTextInput(
                session = session,
                inputId = 'm_est_t',
                value = input$m_est
            )
            
        }
    })
    
    t_est <- reactive({
        if (input$btn_test %% 2 == 0) 
        {
            if(input$radio_test == "One-way unbalanced ANOVA")
            {
                hide("delta1_test")
                hide("delta1_test_t")
                hide("delta2_test")
                hide("delta2_test_t")
            }
            
            hide("sigma_test")
            hide("sigma_test_t")
            hide("FWER_test")
            hide("FWER_test_t")
            hide("nTest_test")
            hide("nTest_test_t")
        }
        if (input$radio_test == "One-way unbalanced ANOVA")
        {
            tab <- data.frame(c(input$radio_test, input$power_test_t, 
                                input$n_test_t, input$maf_test_t, 
                                paste0(input$delta1_test_t, ", ", input$delta2_test_t), 
                                input$sigma_test_t, input$FWER_test_t, input$nTest_test_t))
            rownames(tab)<-c("<strong>Model used</strong>", "<strong>Power level</strong>", 
                             "<strong>Number of subjects needed *</strong>",
                             "<strong>Minor Allele Frequency (MAF)</strong>", 
                             "<strong>Mean differences of gene expression (δ1;δ2)</strong>",
                             "<strong>Standard deviation of gene expression (σ<sub>y</sub>)</strong>",
                             "<strong>Family-wise type I error rate (FWER)</strong>",
                             "<strong>Total number of tests (nTests)</strong>")
        }
        else
        {
            tab <- data.frame(c(input$radio_test, input$power_test_t,
                                input$n_test_t, input$maf_test_t,
                                 input$slope_test_t, input$sigma_test_t, input$FWER_test_t, input$nTest_test_t))
            rownames(tab)<-c("<strong>Model used</strong>", "<strong>Power level</strong>", 
                             "<strong>Number of subjects needed</strong>",
                             "<strong>Minor Allele Frequency (MAF)</strong>", 
                             "<strong>Slope of the regression line (β1)</strong>",
                             "<strong>Standard deviation of gene expression (σ<sub>y</sub>)</strong>",
                             "<strong>Family-wise type I error rate (FWER)</strong>",
                             "<strong>Total number of tests (nTests)</strong>")

        }
        tab
    })
    sc_est <- reactive({
        if (input$btn_scest %% 2 == 0)
        {
            hide("rho_est")
            hide("rho_est_t")
            hide("m_est")
            hide("m_est_t")
            hide("sigma_est")
            hide("sigma_est_t")
            hide("FWER_est")
            hide("FWER_est_t")
            hide("nTest_est")
            hide("nTest_est_t")
        }
        tab <- data.frame(c("Simple linear mixed effects model", input$power_est, ceiling(input$n_est),
                            input$maf_est, input$slope_est,
                            input$sigma_est, input$rho_est, input$m_est,
                            input$FWER_est, input$nTest_est))
        rownames(tab)<-c("<strong>Model used</strong>", "<strong>Power level</strong>", 
                         "<strong>Number of subjects needed *</strong>",
                         "<strong>Minor Allele Frequency (MAF)</strong>", 
                         "<strong>Slope of the regression line (β1)</strong>",
                         "<strong>Standard deviation of gene expression (σ<sub>y</sub>)</strong>",
                         "<strong>Intra-class correlation (ρ)</strong>",
                         "<strong>Number of cells from each subject (m)</strong>",
                         "<strong>Family-wise type I error rate (FWER)</strong>",
                         "<strong>Total number of tests (nTests)</strong>")
        tab
    })
    output$approx <- renderTable(t_est(), rownames = TRUE, colnames = FALSE, sanitize.text.function=function(x){x})
    output$approx1 <- renderTable(sc_est(), rownames = TRUE, colnames = FALSE, sanitize.text.function=function(x){x})
    
    ##Explanation for Sample Estimation
    output$title <- renderUI(tags$div(tags$br(), tags$h4("Tissue eQTL")))
    output$title1 <- renderUI(tags$div(tags$br(), tags$h4("Single-cell eQTL")))
    output$Explanation <- renderUI({
        tags$div(
            tags$br(),
            paste0("* - Minimum number of subjects needed for designated power level."),
            tags$br(),
            tags$br(),
            paste0("More details can be found "),
            HTML("<a onclick=","customHref('about'), target='_blank'>here</a>"),
            paste0("and in the online manual of the function 'powerEQTL.scRNAseq' in R package "),
            tags$a( "‘powerEQTL’", href = "https://CRAN.R-project.org/package=powerEQTL", target="_blank"),
            tags$br(),
            tags$br(),

            )
    })
    
}