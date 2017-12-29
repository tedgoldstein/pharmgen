#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(stringr)
library(shiny)
library(DT)
library(rDGIdb)
library(dplyr)
library(tidyr)


ui <- fluidPage(
    title = "Examples of DataTables",
    sidebarLayout(
        sidebarPanel(
            "Coming soon"

        ),
        mainPanel(
            textAreaInput("Genes", "Genes", value = "TP53 EGFR DDR2", width = "100%", height=30, placeholder = "TP53"),
            tabsetPanel(
                id = 'dataset',
                tabPanel("result",
                         h2("Evidence-Number of articles"),
                         plotOutput("plot1",
                             click = "plot_click",
                             dblclick = "plot_dblclick",
                             hover = "plot_hover",
                             brush = "plot_brush"),
                         verbatimTextOutput("info")),

                tabPanel("table",
                         DT::dataTableOutput("mytable1"))
            )
        )
    )
)


server <- function(input, output) {
    DrugNames = list()

    genes <- reactive({
        str_extract_all(input$Genes, regex("[A-Z0-9a-z-]+", TRUE))[[1]]
    })

    dgidb <- reactive({
        d = NULL
        g = genes()
        if (length(g) > 0) {
            d = queryDGIdb(g)
            d = resultSummary(d)
            d$Drug = toupper(d$Drug)
            d = d %>% group_by(Drug)
            d = d %>% filter(n() > 1)
        }
        d
    })

    output$mytable1 <- DT::renderDataTable({

        DT::datatable(dgidb())
    })

    output$plot1 <- renderPlot({
        plotEvidence(dgidb())
    })

    mergeCommaSets = function(args) {
        x = paste0(args, collapse=",")
        # browser()
        x
    }

    output$info <- renderText({
        d = dgidb()
        dd = d %>%
            group_by(Drug) %>%
            mutate(PMID = mergeCommaSets(PMID))
        ls.str(dd)
    })


    Ndrugs = 10

    plotEvidence <- function(data, ...) {
        n = grep("PMID", colnames(data)) -1
        df = data.frame(Drug=data$Drug, Gene=data$Gene,
                        n = rowSums(data[3:n]))
        df = spread(df, "Drug", "n", fill=0)

        rownames(df) = df$Gene
        df$Gene = NULL
        df = df[, order(colSums(df), decreasing=TRUE)]
        df = df[,1:Ndrugs]
        df = rev(df)

        op <- par(mar = c(5,8,3,1))
        col = rainbow(Ndrugs)
        barplot(as.matrix(df), las = 2, las = 1, horiz = TRUE,
                col=col,
                 ... = ...)
        DrugNames <<- colnames(df)
        legend("bottomright", fill=col, legend=rownames(df))

        invisible()
    }
}

shinyApp(ui, server)
