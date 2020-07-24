library(shiny)
library(shinyjs)
library(shinyBS)
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
    sizeVec.data <- reactive({as.numeric(unlist(strsplit(as.character(input$myntot), ",")))})
    sizeLen.data <- reactive(length(sizeVec.data()))
    # values for sample sizes (single-cell eqtl)
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
    
    ## calculate effect sizes
    effectSize_sc <- reactive({
        ## LME - linear mixed effect model
        sigma2.x = 2*mafVec()[ceiling(mafLen()/2)]*(1-mafVec()[ceiling(mafLen()/2)])
        sigma2.y = input$sigma^2
        
        sigma = sqrt(sigma2.y - input$delta^2*sigma2.x)
        tau = sigma/sqrt(sigma2.x)
        
        round(input$delta/tau,2)
    })
    effectSize_t <- reactive({
        ## SLR - simple linear regression
        sigma2.x = 2*mafVec1()[ceiling(mafLen1()/2)]*(1-mafVec1()[ceiling(mafLen1()/2)])
        sigma2.y = input$sigma1^2
        
        sigma = sqrt(sigma2.y - input$delta1^2*sigma2.x)
        tau = sigma/sqrt(sigma2.x)
        
        round(input$delta1/tau,2)
    })
    
    
    ##author names
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
            
            
            tags$h4("One-way unbalanced ANOVA model"),
            
            # img(src='eQTL_demo_ANOVA.png', align = "center", width = "550px", height = "170px"),
            # tags$br(),
            paste0("Power calculation for eQTL analysis that tests if a SNP is associated to a gene expression level by using unbalanced one-way ANOVA, as the 
                   GTEx consortium used (GTEx consortium, Nature Genetics, 2013)."),
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
            HTML(paste("When δ",tags$sub("1")," = δ",tags$sub("2"), " = δ, we have")),
            helpText(" $$\\lambda = 2 p q  N (\\frac{\\delta}{\\sigma})^2$$"),
            paste0("where p is the minor allele frequency, δ is the mean difference of gene expression levels between genotype classes, and σ is 
                   the standard deviation of gene expression levels in one group of subjects."),
            tags$br(),
            
            # helpText("$$\\text{In this case, the standardized effect size as } \\frac{\\delta}{\\sigma} .$$"),
            tags$br(),

            tags$h4("Simple Linear Regression"),
            
            # img(src='eQTL_demo_SLR.png', align = "center", width = "550px", height = "170px"),
            # tags$br(),
            paste("To test if a SNP is linearly associated with a gene expression level, we use the simple linear regression:"),
            helpText(" $$y_{i} = \\beta_0 + \\beta_1  x_i + \\epsilon_{i}$$"),
            helpText(" $$\\epsilon_{ij} \\text{ ~ } N(0, \\sigma^2_{\\epsilon})$$"),
            helpText("$$i = 1 ... n$$"),
            HTML(paste0("where y", tags$sub("i"), " is the expression level of the gene for the subject i and x", tags$sub("i"), " is the genotype of the i-th 
                    subject by using additive coding (e.g. 0 for AA, 1 for AC, and 2 for CC; A is major allele and C is minor allele). ")),
            # HTML(paste0("β", tags$sub("1")," is the slope of the regression line and can be 
            #         estimated with the following distribution under the alternative hypothesis:")),
            # tags$br(),
            # helpText(" $$\\text{the estimate of the slope } \\hat{\\beta_1} \\text{ ~ } N(\\beta_1, \\frac{\\sigma^2}{( (n-1) \\sigma_x^2)})
            #          \\text{ and the estimate of the intercept } \\hat{\\beta_0} = \\bar{y} - \\hat{\\beta_1} \\bar{x}$$"),
            # helpText("$$\\text{where } \\bar{y} \\text{ is the mean of gene expression, } \\bar{x} \\text{ is the mean of genotype, }$$"),
            # helpText("$$\\sigma_{\\epsilon}^2 = \\sigma^2 - \\beta_1^2 \\sigma_x^2, \\sigma^2 \\text{ is the variance of the outcome, } \\sigma_x^2 
            #          \\text{ is the variance of predictor, and } \\sigma_{\\epsilon}^2 \\text{ is the variance of the random error.}$$"),
            
            paste0("We can derive the power calculation formula as: "),
            helpText("$$power = 1 - T_{n-2, \\lambda}[t_{n-2}(\\alpha/2)] + T_{n-2, \\lambda}[-t_{n-2}(\\alpha/2)]$$"),
            HTML(paste0("where T", tags$sub("n−2,λ"), "(a) is the value at a of the cumulative distribution function of non-central t distribution with (n − 2)
                    degrees of freedom and non-centrality parameter λ. λ can be written as: ")),
            helpText("$$ \\lambda = \\frac{\\beta_1}{\\sqrt{(\\sigma_y^2 - \\beta_1^2 2pq) / [(n-1)2pq]}} $$"),
            paste0("p is the minor allele frequency (MAF) and q = 1-p."),
            # helpText("$$\\text{The unstandardized effect size as } \\beta_1", "\\text{, and
            #                 the standardized effect size equal to } \\frac{\\beta_1}{\\tau} \\text{, where }\\tau = \\frac{\\sigma}{\\sigma_x} \\text{ is the
            #                 standard deviation of random error}$$"),
            # helpText("$$\\text{ normalized by the standard deviation of predictor}.$$"),
            
            tags$br(),
            
            tags$h4("Simple linear mixed effects model"),
            
            paste0("For single-cell eQTL, we assume the following simple linear mixed effects model for each (SNP, gene) pair to 
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
            helpText("$$\\text{where } a=\\frac{\\sqrt{2p(1-p)}}{\\sigma_y} \\text{ and } b = \\frac{\\beta_1 \\sqrt{m(n-1)}}{\\sqrt{1+(m-1)\\rho}}$$"),
            HTML(paste0("z",tags$sub("α∗/2"), " is the upper 100α∗/2 percentile of the standard normal distribution,
                        α∗ = α/nTests, nTests is the number of (SNP, gene) pairs, ρ is the intraclass correlation,.")),
            
            # helpText("$$\\text{We define the unstandardized effect size as } \\beta_1", "\\text{, and
            #                 the standardized effect size equal to } \\frac{\\beta_1}{\\tau} \\text{, where }\\tau = \\frac{\\sigma_{\\epsilon}}{\\sigma_x} \\text{ is the
            #                 standard deviation of random error}$$"),
            # helpText("$$\\text{ normalized by the standard deviation of predictor}.$$"),
            tags$br(),
            tags$br(),
            
            tags$b("References"),
            tags$br(),
            
            paste0("O'Brien RG Muller KE. Unified power analysis for t-tests through multivariate hypotheses. Applied Analysis of Variance in the Behavioral Sciences, 297–344, 1993. "),
            tags$br(),
            paste0("The GTEx Consortium. The Genotype-Tissue Expression (GTEx) project. Nature Genetics, 45:580-585, 2013."),
            tags$br(),
            tags$br(),
            
            tags$h4("Developer and Bug Report"),
            paste0("This R shiny application is developed by Xiaoqi Li (xli778@wisc.edu), under the supervisions of Drs. 
                    Xianjun Dong and Weiliang Qiu."),
            tags$br(),
            
            paste0("Thank you very much for using our power calculator for bulk tissue 
                   and single-cell eQTL analysis! Please let us know what do you think."),
            tags$a("Click here!", href = "https://gitreports.com/issue/bwh-bioinformatics-hub/Rshiny_PowerAnalysis"),
            tags$br(),
            tags$br(),
            
            tags$h4("To cite our work, please use citation below:"),
            paste0("Dong X, Li X, Chang T, Weiss S, and Qiu W. powerEQTL: an R package and R shiny application for calculating 
                    sample size and power of bulk tissue and single-cell eQTL analysis. manuscript. (2020)"),
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
        colors <- hcl.colors(sizeLen.data(),"Set 2")
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
            colors <- hcl.colors(sizeLen.data(),"Set 2")
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
            show("delta2")
            output$description2 <- renderUI({
                tags$div(
                    withMathJax(),
                    
                    # tags$br(),
                    
                    paste0("The figure above shows the statistical power of an eQTL analysis as a function of the minor 
                            allele frequency (MAF) with different sample size. Unbalanced one-way ANOVA test was used 
                            to test the potential nonlinear association of genotype classes to expression level, as Ref.
                            (GTEx consortium, Nature Genetics, 2013). We model the tansformed expression data (e.g. using
                            log transformation) as normally distributed with a standard deviation of ", input$sigma1, " (we assume
                            that each genotype class of subjects have the same standard deviation). This level of noise 
                            can be based on estimates from previous study or pilot data. "),
                    # paste0("The effect size depends both on the MAF of the SNPs and the actual log expression change between
                    #         genotype classes (denoted by δ1 and δ2). In this case, when δ1 = δ2 = δ, δ = ", input$delta1, " and 
                    #         standard deviation = ", input$sigma1, " is equivalent to detecting a ", input$delta1/input$sigma1,
                    #         "-fold (", input$delta1/input$sigma1, " = ",  input$delta1, "/", input$sigma1, ") of expression 
                    #         level chance across the genotype classes. "),
                    # paste0(" The standardized effect size is ", 
                    #         input$delta1/input$sigma1," in this example."),
                    tags$br(),
                    tags$br(),
                    
                    paste0("This plot uses the function of ‘powerEQTL.ANOVA’ in R package "),
                    tags$a( "‘powerEQTL’", href = "https://CRAN.R-project.org/package=powerEQTL"),
                    paste0(". More details can be found "),
                    HTML("<a onclick=","customHref('about')>here</a>"),
                    paste0("and in the online manual of the function ‘powerEQTL.ANOVA’ in R package "),
                    tags$a( "‘powerEQTL’", href = "https://CRAN.R-project.org/package=powerEQTL"),
                    
                    tags$br(),
                )
            })
        }
        else
        {
            updateNumericInput(session, inputId = "delta1", label = "Slope of the regression line (β1)", value = input$delta1)
            hide("delta2")
            output$description2 <- renderUI({
                tags$div(
                    
                    withMathJax(),
                    
                    # tags$br(),
                    
                    paste0("The figure above shows the statistical power of an eQTL analysis as a function of the minor
                            allele frequency (MAF) with different sample size.  Simple linear regression model was used 
                            to test the linear association between gene expression and genotype.
                            In this example, we model the expression data as log-normally distributed with a log standard
                            deviation of ", input$sigma1, ". This level of noise can be based on estimates from previous
                            study or pilot data. "),
                    # paste0(" In this example, the unstandardized effect size as β1 = ", input$delta1,
                    #         ", and the standardized effect size equal to ", effectSize_t(),"."),
                    tags$br(),
                    tags$br(),
                    
                    paste0("This plot uses the function of ‘powerEQTL.SLR’ in R package "),
                    tags$a( "‘powerEQTL’", href = "https://CRAN.R-project.org/package=powerEQTL"),
                    paste0(". More details can be found "),
                    HTML("<a onclick=","customHref('about')>here</a>"),
                    paste0("and in the online manual of the function ‘powerEQTL.SLR’ in R package "),
                    tags$a( "‘powerEQTL’", href = "https://CRAN.R-project.org/package=powerEQTL"),
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
            
            legend(input$pos, title="Sample size (n)", legend = sizeVec.data(),
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
                    frequency (MAF) with different sample size.  Simple linear mixed effect regression model was used 
                    to characterize the linear association between gene expression and genotype. We assume expression 
                    levels are normally distributed after appropriate pre-processing. In this example, We model the 
                    expression data as log-normally distributed with a log standard deviation of ", input$sigma, " 
                    and an intra-class correlation of ", input$rho, ". This level of noise can be based on estimates 
                    from previous study or pilot data. ")),
                    
        # paste0("The standardized effect size depends on the MAF of the SNPs, the standard deviation of gene expression, and the slope 
        # of linear regression (denoted by β).In this case, the unstandardized effect size is β =", input$delta, ", and the standardized effect size 
        #             equal to ", effectSize_sc(),"."),
        
        tags$br(),
        tags$br(),  
        
        paste0("This plot uses function ‘powerEQTL.scRNAseq’ in R package "),
        tags$a( "‘powerEQTL’", href = "https://CRAN.R-project.org/package=powerEQTL"),
        paste0(". More details can be found "),
        HTML("<a onclick=","customHref('about')>here</a>"),
        paste0("and in the online manual of the function ‘powerEQTL.scRNAseq’ in R package "),
        tags$a( "‘powerEQTL’", href = "https://CRAN.R-project.org/package=powerEQTL"),
        tags$br(),
         
        
        # tags$b("References"),
        # tags$br(),
        # 
        # paste0("Dong X, Li X, Chang T, Weiss S, and Qiu W. powerEQTL: an R package and R shiny application for calculating 
        #        sample size and power of bulk tissue and single-cell eQTL analysis. manuscript. (2020)")
        )
    })
    
    ## Power calculator for single-cell eQTL
    ## hidden input options
    observeEvent(input$btn_scest, {
        toggle("sigma_est")
        toggle("FWER_est")
        toggle("nTest_est")
        toggle("rho_est")
        toggle("m_est")
        if (input$btn_scest %% 2 == 1) {
            txt <- "Less Options"
        } else {
            txt <- "More Options"
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
    })
    
    ## Power calculator for tissue eQTL
    # Reactive Input Title for delta for different models
    observeEvent(input$radio_test,{
        if(input$radio_test == "One-way unbalanced ANOVA")
        {
            hide("slope_test")
            if (input$btn_test %% 2 == 1) 
            {
                toggle("delta1_test")
                toggle("delta2_test")
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
            
            updateSliderInput(session, inputId = "n_test", label = "Number of subjects needed",
                              value = as.numeric(n), max = ceiling(as.numeric(max)))
            
            output$Explanation1 <- renderUI({
                tags$div(
                    tags$br(),
                    paste0("The table uses function 'powerEQTL.ANOVA' in R package "),
                    tags$a( "‘powerEQTL’", href = "https://CRAN.R-project.org/package=powerEQTL"),
                    paste0(". More details can be found "),
                    HTML("<a onclick=","customHref('about')>here</a>"),
                    paste0("and in the online manual of the function 'powerEQTL.ANOVA' in R package "),
                    tags$a( "‘powerEQTL’", href = "https://CRAN.R-project.org/package=powerEQTL"),
                    tags$br(),
                    tags$br(),
                )
            })
        }
        else
        {
            show("slope_test")
            hide("delta1_test")
            hide("delta2_test")
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
            updateSliderInput(session, inputId = "n_test", label = "Number of subjects needed",
                              value = as.numeric(n), max = ceiling(as.numeric(max)))
            
            output$Explanation1 <- renderUI({
                tags$div(
                    tags$br(),
                    paste0("The table uses function 'powerEQTL.SLR' in R package "),
                    tags$a( "‘powerEQTL’", href = "https://CRAN.R-project.org/package=powerEQTL"),
                    paste0(". More details can be found "),
                    HTML("<a onclick=","customHref('about')>here</a>"),
                    paste0("and in the online manual of the function 'powerEQTL.SLR' in R package "),
                    tags$a( "‘powerEQTL’", href = "https://CRAN.R-project.org/package=powerEQTL"),
                    tags$br(),
                    tags$br(),
                    
                    # tags$b("References"),
                    # tags$br(),
                    # 
                    # paste0("Dong X, Li X, Chang T, Weiss S, and Qiu W. powerEQTL: an R package and R shiny application for calculating 
                    #    sample size and power of bulk tissue and single-cell eQTL analysis. manuscript. (2020)"),
                    # tags$br(),
                    # paste0("Dupont, W.D. and Plummer, W.D.. Power and Sample Size Calculations for Studies Involving Linear Regression. Controlled Clinical Trials. 1998;19:589-601."),
                    # tags$br(),
                    # paste0("Lonsdale J and Thomas J, et al. The Genotype-Tissue Expression (GTEx) project. Nature Genetics, 45:580-585, 2013.")
                    
                )
            })
            
            
        }

    })
    
    ## hidden input options
    observeEvent(input$btn_test, {
        if(input$radio_test == "One-way unbalanced ANOVA")
        {
            toggle("delta1_test")
            toggle("delta2_test")
        }
       
        toggle("sigma_test")
        toggle("FWER_test")
        toggle("nTest_test")
        if (input$btn_test %% 2 == 1) {
            txt <- "Less Options"
        } else {
            txt <- "More Options"
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
            
            updateSliderInput(session, inputId = "n_test", label = "Number of subjects needed",
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
            updateSliderInput(session, inputId = "n_test", label = "Number of subjects needed",
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
            
            updateSliderInput(session, inputId = "n_test", label = "Number of subjects needed",
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
            updateSliderInput(session, inputId = "n_test", label = "Number of subjects needed",
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
            updateSliderInput(session, inputId = "n_test", label = "Number of subjects needed",
                              value = ceiling(as.numeric(n)), max = ceiling(as.numeric(max)))
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
            
            updateSliderInput(session, inputId = "n_test", label = "Number of subjects needed",
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
            updateSliderInput(session, inputId = "n_test", label = "Number of subjects needed",
                              value = ceiling(as.numeric(n)), max = ceiling(as.numeric(max)))
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
            
            updateSliderInput(session, inputId = "n_test", label = "Number of subjects needed",
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
            updateSliderInput(session, inputId = "n_test", label = "Number of subjects needed",
                              value = ceiling(as.numeric(n)), max = ceiling(as.numeric(max)))
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
            
            updateSliderInput(session, inputId = "n_test", label = "Number of subjects needed",
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
            updateSliderInput(session, inputId = "n_test", label = "Number of subjects needed",
                              value = ceiling(as.numeric(n)), max = ceiling(as.numeric(max)))
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
            
            updateSliderInput(session, inputId = "n_test", label = "Number of subjects needed",
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
            updateSliderInput(session, inputId = "n_test", label = "Number of subjects needed",
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
            
            updateSliderInput(session, inputId = "n_test", label = "Number of subjects needed",
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
            updateSliderInput(session, inputId = "n_test", label = "Number of subjects needed",
                              value = ceiling(as.numeric(n)), max = ceiling(as.numeric(max)))
        }
    })
    
    t_est <- reactive({
        if (input$radio_test == "One-way unbalanced ANOVA")
        {
            tab <- data.frame(c(input$radio_test, input$power_test, 
                                input$n_test, input$maf_test, 
                                paste0(input$delta1_test, ", ", input$delta2_test), 
                                input$sigma_test, input$FWER_est, input$nTest_est))
            rownames(tab)<-c("<strong>Model used</strong>", "<strong>Power level</strong>", 
                             "<strong>Number of subjects needed</strong>",
                             "<strong>Minor Allele Frequency (MAF)</strong>", 
                             "<strong>Mean differences of gene expression (δ1;δ2)</strong>",
                             "<strong>Standard deviation of gene expression (σ<sub>y</sub>)</strong>",
                             "<strong>Family-wise type I error rate (FWER)</strong>",
                             "<strong>Total number of tests (nTests)</strong>")
        }
        else
        {
            tab <- data.frame(c(input$radio_test, input$power_test,
                                input$n_test, input$maf_test,
                                 input$slope_test, input$sigma_test, input$FWER_est, input$nTest_est))
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
        tab <- data.frame(c("Simple linear mixed effects model", input$power_est, ceiling(input$n_est),
                            input$maf_est, input$slope_est,
                            input$sigma_est, input$rho_est, input$m_est,
                            input$FWER_est, input$nTest_est))
        rownames(tab)<-c("<strong>Model used</strong>", "<strong>Power level</strong>", 
                         "<strong>Number of subjects needed</strong>",
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
            paste0("The table uses function 'powerEQTL.scRNAseq' in R package "),
            tags$a( "‘powerEQTL’", href = "https://CRAN.R-project.org/package=powerEQTL"),
            paste0(". More details can be found "),
            HTML("<a onclick=","customHref('about')>here</a>"),
            paste0("and in the online manual of the function 'powerEQTL.scRNAseq' in R package "),
            tags$a( "‘powerEQTL’", href = "https://CRAN.R-project.org/package=powerEQTL"),
            tags$br(),
            tags$br(),
            
            # tags$b("References"),
            # tags$br(),
            # 
            # paste0("Dong X, Li X, Chang T, Weiss S, and Qiu W. powerEQTL: an R package and R shiny application for calculating 
            #    sample size and power of bulk tissue and single-cell eQTL analysis. manuscript. (2020)"),
            # tags$br(),
            # paste0("Dupont, W.D. and Plummer, W.D.. Power and Sample Size Calculations for Studies Involving Linear Regression. Controlled Clinical Trials. 1998;19:589-601."),
            # tags$br(),
            # paste0("Lonsdale J and Thomas J, et al. The Genotype-Tissue Expression (GTEx) project. Nature Genetics, 45:580-585, 2013.")
            
            )
    })
    
}