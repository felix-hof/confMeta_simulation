library(shiny)
library(ggplot2)

## UI
ui <- fluidPage(
    titlePanel("Simulation study from Held, Hofmann, XXXX (2023)"),
    tabsetPanel(
        tabPanel("Visualize results",
                 br(),
                 sidebarLayout(
                     sidebarPanel(
                         selectInput(inputId = "outcome", label = "Performance measure",
                                     choices = c("Coverage (true overall effect)" = "coverage_true",
                                                 "Coverage (true study effects - proportion of effects)" = "coverage_effects",
                                                 "Coverage (true study effects - all included)" = "coverage_effects_all",
                                                 "Coverage (true study effects - at least one included)" = "coverage_effects_min1",
                                                 "Coverage (prediction of new study effect)" = "coverage_prediction",
                                                 "Number of intervals" = "n",
                                                 "Interval score" = "score",
                                                 "Interval width" = "width")
                                     ),
                         br(),
                         h3("Data generating process"),
                         selectInput(inputId = "bias", label = "Bias",
                                     choices = c("none",
                                                 "moderate",
                                                 "strong")
                                     ),
                         selectInput(inputId = "dist", label = "Distribution",
                                     choices = c("Gaussian",
                                                 "t")
                                     ),
                         selectInput(inputId = "model", label = "Simulation model",
                                     choices = c("additive",
                                                 "multiplicative")
                                     ),
                         selectInput(inputId = "nlarge", label = "Number of large studies",
                                     choices = c(0, 1, 2)
                                     ),
                         br(),
                         h3("Plotting options"),
                         p("TODO: Add options")
                     ),
                     mainPanel(
                         plotOutput(outputId = "simulationplot")
                     )
                 )
                 ),
        tabPanel("Simulation protocol",
                 br(),
                 p("lol")
                 )
    )
)

## server
server <- function(input, output) {
    output$simulationplot <- renderPlot({
        ## TODO add code for creating plots
        ggplot(iris, aes(x = Sepal.Width, y = Petal.Length)) +
            geom_point()
    })
}

## Run the app
shinyApp(ui = ui, server = server)
