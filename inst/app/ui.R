#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

    # Application title
    titlePanel("Gene Abundance by Tissue and Study Arm"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            textInput(inputId = "gene", label = "Gene" , value = "Cd3d"),
            out = h5("Cd3 = Cd3d; Cd8 = Cd8a; Cd4 = Cd4"),
            actionButton(inputId = "find_gene", label = "Run")
        ),
        
        mainPanel(
        tabsetPanel(
            tabPanel("Plot", plotOutput("plotOut")),
            tabPanel("Log2TPM Data", downloadButton('downloadData', 'Download CSV'), DT::dataTableOutput("genetable")),
            tabPanel("Treatment Data", DT::dataTableOutput("trtable"))
        )
    )
    )
))
