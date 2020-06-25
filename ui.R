library(shiny)

basicPage(
    titlePanel("Power Analysis"),
    # fluidRow(column(width = 12, htmlOutput("description"))),
    tabsetPanel(type = "tabs",
                tabPanel("single-cell eQTL", 
                         sidebarLayout(
                             sidebarPanel(
                                 
                                 textInput(inputId = "subjects1", label = "Subject names (separated by ',')",
                                           value = "DA neurons, Glu neurons, GABA neuron, Astrocytes, Oligodendrocytes, Microglia"),
                                 textInput(inputId = "myntot1", label = "Sample sizes (separated by ',')", 
                                           value = "50, 100, 150, 200, 250, 300"),
                                 numericInput(inputId = "m", label = "Number of cells from each subject (m)", 
                                              value = 640, min = 0),
                                 numericInput(inputId = "sigma", label = "Standard deviation of gene expression levels in one group of subjects.", 
                                              value = 0.29),
                                 numericInput(inputId = "rho", label = "Intra-class correlation (range [0, 1])", 
                                              value = 0.5, min = 0, max = 1),
                                 numericInput(inputId = "alpha", label = "Family-wise type I error rate", 
                                              value = 0.05),
                                 numericInput(inputId = "nTest", label = "Total number of tests", 
                                              value = 10e5),
                                 numericInput(inputId = "delta", label = "Mean difference of gene expression levels between groups (slope)", 
                                              value = 0.29*1.5),

                                 submitButton("submit"),
                                 
                                 htmlOutput("description1")
                             ),
                             mainPanel(
                                 
                                 plotOutput(outputId = "cell", click = "plot_click", width = 800, height = 500),

                                 downloadButton('export', "Download Report"),
                                 img(src='eqtl.png', align = "middle", width = 480, height = 200),
                                 htmlOutput("description"),
                                 
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
                                 numericInput(inputId = "delta1", label = "Mean difference of gene expression levels between groups (slope)", 
                                              value = 0.29*1.5),
                                 numericInput(inputId = "alpha1", label = "Family-wise type I error rate", 
                                              value = 0.05),
                                 numericInput(inputId = "nTest1", label = "Total number of tests", 
                                              value = 10e5),
                                 radioButtons(inputId = "radio", label = "Model",
                                              choices = c("One-way unbalanced anova", "Simple linear regression")),
                                 submitButton("submit"),

                                 htmlOutput("description3")
                             ),
                             mainPanel(
                                 
                                 plotOutput(outputId = "tissue", click = "plot_click1", width = 800, height = 500),
                                 downloadButton('export2', "Download Report"),
                                 
                                 img(src='eqtl.png', align = "middle", width = 480, height = 200),
                                 htmlOutput("description2"),
                                 
                                 
                                 
                             ))),
                tabPanel("Sample Size Estimation",
                         sidebarLayout(
                             sidebarPanel(
                                 
                                 numericInput(inputId = "maf2", label = "Minor allele frequencies (between 0 and 0.5)",
                                              value = 0.02),
                                 numericInput(inputId = "power", label = "Desired Power level",
                                              value = 0.8, min = 0, max = 1),
                                 numericInput(inputId = "eSize", label = "eQTL effect size", 
                                              value = 1),
                                 numericInput(inputId = "rho3", label = "Intra-class correlation (range [0, 1]; single-cell)", 
                                              value = 0.5, min = 0, max = 1),
                                 numericInput(inputId = "m3", label = "Number of cells from each subject (m; single-cell)", 
                                              value = 640, min = 0),
                                 
                                 
                                 submitButton("submit")
                                ),
                             mainPanel(
                                 htmlOutput("title"),
                                 tableOutput("approx"),
                                 htmlOutput("title1"),
                                 tableOutput("approx1"),
                                 htmlOutput("Explanation")
                                 )
                             )
                         )
                # tabPanel("Single-cell eQTL - Sample Size Estimation",
                #          sidebarLayout(
                #              sidebarPanel(
                #                  
                #                  numericInput(inputId = "maf3", label = "Minor allele frequencies (between 0 and 0.5)",
                #                               value = 0.02),
                #                  numericInput(inputId = "power1", label = "Desired Power level",
                #                               value = 0.8, min = 0, max = 1),
                #                  numericInput(inputId = "delta3", label = "eQTL effect size", 
                #                               value = 0.29*1.5),
                #                  numericInput(inputId = "rho3", label = "Intra-class correlation (range [0, 1])", 
                #                               value = 0.5, min = 0, max = 1),
                #                  numericInput(inputId = "m3", label = "Number of cells from each subject (m)", 
                #                               value = 640, min = 0),
                #                  
                #                  submitButton("submit")
                #              ),
                #              mainPanel(
                #                  tableOutput("approx1")
                #              )
                #          )
                # )
                
                #htmlOutput("summary"),
    )
)

