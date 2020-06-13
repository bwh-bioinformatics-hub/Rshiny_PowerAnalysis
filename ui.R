library(shiny)

basicPage(
    titlePanel("Power Analysis"),
    # fluidRow(column(width = 12, htmlOutput("description"))),
    tabsetPanel(type = "tabs",
                tabPanel("single-cell eQTL", 
                         sidebarLayout(
                                sidebarPanel(
                                        numericInput(inputId = "n", label = "Number of subjects (n)", 
                                                     value = 100, min = 0),
                                        numericInput(inputId = "m", label = "Number of cells from each subject (m)", 
                                                     value = 640, min = 0),
                                        numericInput(inputId = "sigma", label = "Standard deviation of gene expression levels in one group of subjects.", 
                                                     value = 0.29),
                                        numericInput(inputId = "rho", label = "Intra-class correlation (range [0, 1])", 
                                                     value = 0.5, min = 0, max = 1),
                                        numericInput(inputId = "alpha", label = "Family-wise type I error rate", 
                                                     value = 5.4e-8, min = 0, max = 1),
                                        textInput(inputId = "MAF", label = "Minor allele frequencies (between 0 and 0.5; separated by ',')", 
                                                     value = "0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1"),
                                        submitButton("submit"),
                                        plotOutput(outputId = "effectSize", width = 350, height = 200),
                                        htmlOutput("description1")
                                    ),
                                mainPanel(
                                    
                                    plotOutput(outputId = "cell", click = "plot_click", width = 800, height = 500),
                                    # verbatimTextOutput("cell_info"),
                                    downloadButton('export', "Download Report"),
                                    htmlOutput("description")
                                    )
                                )
                ),
                tabPanel("tissue eQTL", 
                         sidebarLayout(
                                sidebarPanel(
                                    
                                        textInput(inputId = "subjects", label = "Subject names (separated by ',')",
                                                  value = "DA neurons, Glu neurons, GABA neuron, Astrocytes, Oligodendrocytes, Microglia"),
                                        textInput(inputId = "myntot", label = "Sample sizes (separated by ',')", 
                                                     value = "50, 100, 150, 200, 250, 300"),
                                        numericInput(inputId = "sigma1", label = "Standard deviation of gene expression levels in one group of subjects.", 
                                                     value = 0.29),
                                        numericInput(inputId = "delta", label = "eQTL effect size", 
                                                     value = 0.29*1.5),
                                        numericInput(inputId = "alpha1", label = "Family-wise type I error rate = 0.05/nTests = 0.05/(number of genes * number of SNPs)", 
                                                     value = 5.4e-8, min = 0, max = 1),
                                        submitButton("submit"),
                                        plotOutput(outputId = "effectSize1", width = 350, height = 200),
                                        htmlOutput("description3")
                                ),
                                mainPanel(
                                    
                                        plotOutput(outputId = "tissue", click = "plot_click1", width = 800, height = 500),
                                        # verbatimTextOutput("tissue_info"),
                                        downloadButton('export2', "Download Report"),
                                        htmlOutput("description2")
                                    
                                    
                                )))
                                    
                #htmlOutput("summary"),
        )
    )
