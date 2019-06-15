library(shiny)
library(bSims)

ui <- navbarPage("bSims (H)",
  tabPanel("Initialize",
    column(6,
      plotOutput(outputId = "plot_ini")),
    column(6,
      sliderInput("seed", "Random seed", 0, 100, 0, 1)
    )
  ),
  tabPanel("Populate",
    column(6,
      plotOutput(outputId = "plot_pop")),
    column(6,
      sliderInput("D", "Density", 0, 10, 1, 0.1),
      radioButtons("spfun", "Spatial pattern",
        c("Random"="random", "Regular"="regular",
          "Clustered"="clustered"))
    )
  ),
  tabPanel("Animate",
    column(6,
      plotOutput(outputId = "plot_ani")),
    column(6,
      sliderInput("phi1", "Vocal rate (group 1)", 0, 10, 0.5, 0.1),
      sliderInput("phi2", "Vocal rate (group 2)", 0, 10, 0, 0.1),
      sliderInput("mix", "Mixture (group 1)", 0, 1, 1, 0.05),
      sliderInput("phim", "Movement rate", 0, 10, 1, 0.1),
      sliderInput("SDm", "Movement SD", 0, 1, 0, 0.05)
    )
  ),
  tabPanel("Detect",
    column(6,
      plotOutput(outputId = "plot_det")),
    column(6,
      sliderInput("tau", "Detection parameter", 0, 2, 1, 0.1),
      radioButtons("dfun", "Distance function",
        c("Half Normal"="halfnormal", "Negative Exponential"="negexp")),
      sliderInput("repel", "Repel distance", 0, 2, 0, 0.1)
    )
  ),
  tabPanel("Transcribe",
    column(6,
      plotOutput(outputId = "plot_tra")),
    column(6,
      sliderInput("y", "y", 0, 10, 1, 0.1)
    )
  ),
  tabPanel("Estimate",
    column(6,
      plotOutput(outputId = "plot_est")),
    column(6,
      sliderInput("q", "q", 0, 10, 1, 0.1)
    )
  )
)


# plotOutput(outputId = "distPlot")

server <- function(input, output) {
  l <- reactive({
    set.seed(input$seed)
    bsims_init(extent = 10)
  })
  a <- reactive({
    xy_fun <- switch(input$spfun,
      "random"=function(d) ifelse(d > 0, 1, 1),
      "regular"=function(d) 1-exp(-d^2/2^2),
      "clustered"=function(d)
        pmax(ifelse(d < 0.1, 1, 0), 0.5*(1-exp(-d^2/2^2))))
    margin <- switch(input$spfun,
      "random"=0,
      "regular"=2,
      "clustered"=5)
    bsims_populate(l(),
      density = input$D,
      xy_fun = xy_fun,
      margin = margin)
  })
  b <- reactive({
    bsims_animate(a(),
      duration = 10,
      vocal_rate = c(input$phi1, input$phi2),
      move_rate = input$phim,
      movement = input$SDm,
      mixture = c(input$mix, 1-input$mix))
  })
  d <- reactive({
    dfun <- switch(input$dfun,
      "halfnormal"=function(d, tau) exp(-d^2/tau^2),
      "negexp"    =function(d, tau) exp(-d/tau^2))
    bsims_detect(b(),
      xy = c(0, 0),
      tau = input$tau,
      dist_fun = NULL,
      repel = input$repel,
      vocal_only = TRUE)
  })

  output$plot_ini <- renderPlot({
    op <- par(mar=c(0,0,0,0))
    plot(l())
    par(op)
  })
  output$plot_pop <- renderPlot({
    op <- par(mar=c(0,0,0,0))
    plot(a())
    par(op)
  })
  output$plot_ani <- renderPlot({
    op <- par(mar=c(0,0,0,0))
    plot(b())
    par(op)
  })
  output$plot_det <- renderPlot({
    op <- par(mar=c(0,0,0,0))
    plot(d())
    par(op)
  })
  output$plot_tra <- renderPlot({
    op <- par(mar=c(0,0,0,0))
    plot(d())
    par(op)
  })
  output$plot_est <- renderPlot({
    op <- par(mar=c(0,0,0,0))
    plot(d())
    par(op)
  })
}

shinyApp(ui = ui, server = server)
