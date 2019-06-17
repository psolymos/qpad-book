library(shiny)
library(detect)
library(bSims)
source("../functions.R")

EXTENT <- 10
DURATION <- 10
TINT <- list(
  "0-10 min"=c(10),
  "0-3-5-10 min"=c(3, 5, 10),
  "0-11-2-3 min"=c(1, 2, 3)
)
RINT <- list(
  "0-Inf m"=c(Inf),
  "0-50-100-Inf m"=c(0.5, 1, Inf),
  "0-50-100 m"=c(0.5, 1),
  "0-50 m"=c(0.5)
)

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
      plotOutput(outputId = "plot_pop"),
      plotOutput(outputId = "plot_spfun")
    ),
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
      plotOutput(outputId = "plot_det"),
      plotOutput(outputId = "plot_dfun")
    ),
    column(6,
      sliderInput("tau", "Detection parameter (SD)", 0, 5, 1, 0.25),
      sliderInput("bpar", "Hazard rate parameter (SD)", 0, 5, 1, 0.5),
      radioButtons("dfun", "Distance function",
        c("Half Normal"="halfnormal",
          "Negative Exponential"="negexp",
          "Hazard rate"="hazrate")),
      sliderInput("repel", "Repel distance", 0, 2, 0, 0.1),
      checkboxInput("onlyv", "Vocalizations only", TRUE)
    )
  ),
  tabPanel("Transcribe",
    fluidRow(
      column(6,
        plotOutput(outputId = "plot_tra")
      ),
      column(6,
        selectInput("tint", "Time intervals", names(TINT)),
        selectInput("rint", "Distance intervals", names(RINT)),
        sliderInput("derr", "Distance error", 0, 1, 0, 0.1),
        checkboxInput("only1", "1st detection only", TRUE)
      )
    ),
    fluidRow(
      column(6,
        h4("Removal table"),
        tableOutput(outputId = "table_rem")
      ),
      column(6,
        h4("Multiple-visits"),
        tableOutput(outputId = "table_vis")
      )
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
  dis <- seq(0, 10, 0.01)
  l <- reactive({
    set.seed(input$seed)
    bsims_init(extent = EXTENT)
  })
  xy_fun <- reactive({
    switch(input$spfun,
      "random"=function(d) rep(1, length(d)),
      "regular"=function(d)
        (1-exp(-d^2/1^2) + dlnorm(d, 2)/dlnorm(2,2)) / 2,
      "clustered"=function(d)
        exp(-d^2/1^2) + 0.5*(1-exp(-d^2/4^2))
      )
  })
  a <- reactive({
    margin <- switch(input$spfun,
      "random"=0,
      "regular"=2,
      "clustered"=5)
    bsims_populate(l(),
      density = input$D,
      xy_fun = xy_fun(),
      margin = margin)
  })
  b <- reactive({
    bsims_animate(a(),
      duration = DURATION,
      vocal_rate = c(input$phi1, input$phi2),
      move_rate = input$phim,
      movement = input$SDm,
      mixture = c(input$mix, 1-input$mix))
  })
  dfun <- reactive({
    switch(input$dfun,
      "halfnormal"=function(d, tau) exp(-(d/tau)^2),
      "negexp"    =function(d, tau) exp(-d/tau),
      "hazrate"   =function(d, tau) 1-exp(-(d/tau)^-input$bpar)
    )
  })
  d <- reactive({
    bsims_detect(b(),
      xy = c(0, 0),
      tau = input$tau,
      dist_fun = dfun(),
      repel = input$repel,
      vocal_only = input$onlyv)
  })
  m <- reactive({
    bsims_transcribe(d(),
      tint = TINT[[input$tint]],
      rint = RINT[[input$rint]],
      first_only = input$only1,
      error = input$derr
    )
  })
  e <- reactive({
    Ydur <- matrix(colSums(m()$removal), 1)
    Ddur <- matrix(TINT[[input$tint]], 1)
    Ydis <- matrix(rowSums(m()$removal), 1)
    Ddis <- matrix(RINT[[input$rint]], 1)
    if (length(TINT[[input$tint]]) > 1) {
      Mrem <- cmulti.fit(Ydur, Ddur, type="rem")
      Mmix <- cmulti.fit(Ydur, Ddur, type="mix")
    } else {
      Mrem <- NULL
      Mmix <- NULL
    }
    if (length(RINT[[input$rint]]) > 1) {
      Mdis <- cmulti.fit(Ydis, Ddis, type="dis")
    } else {
      Mdis <- NULL
    }
    list(
      Ydur=Ydur, Ddur=Ddur,
      Ydis=Ydis, Ddis=Ddis,
      Mrem=Mrem, Mmix=Mmix, Mdis=Mdis)
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
  output$plot_spfun <- renderPlot({
    plot(dis, xy_fun()(dis), type="l", col=4,
      ylim=c(0,1), xlab="Distance", ylab="P(acceptence)")
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
  output$plot_dfun <- renderPlot({
    plot(dis, dfun()(dis, input$tau), type="l", col=4,
      ylim=c(0,1), xlab="Distance", ylab="P(detection)")
  })
  output$plot_tra <- renderPlot({
    op <- par(mar=c(0,0,0,0))
    plot(m(), tlim=c(0, max(TINT[[input$tint]])))
    rr <- RINT[[input$rint]]
    rr <- rr[is.finite(rr)]
    if (length(rr) > 0) {
      if (any(is.infinite(RINT[[input$rint]]))) {
        polygon(0.5*EXTENT*c(-1,-1,1,1), 0.5*EXTENT*c(-1,1,1,-1),
          border=NA, col="#ff000033")
      } else {
        draw_ellipse(0, 0, max(rr), max(rr),
          border=NA, col="#ff000033")
      }
      draw_ellipse(rep(0, length(rr)), rep(0, length(rr)), rr, rr,
        border=2)
    } else {
      polygon(0.5*EXTENT*c(-1,-1,1,1), 0.5*EXTENT*c(-1,1,1,-1),
        border=NA, col="#ff000033")
    }
    tt <- TINT[[input$tint]]
    tt <- EXTENT * TINT[[input$tint]] / DURATION
    tt <- c(0, tt) * 0.8 - (EXTENT * 0.4)
    if (max(TINT[[input$tint]]) < DURATION)
      polygon(EXTENT*0.4*c(-1,-1,1,1), EXTENT*c(0.4, 0.45, 0.45, 0.4),
        border=2, col=NA, lty=2)
    polygon(EXTENT*0.4*c(-1,-1,1,1), EXTENT*c(0.4, 0.45, 0.45, 0.4),
      border=2, col=NA, lty=2)
    for (i in 2:length(tt)) {
      polygon(tt[c(i-1, i-1, i, i)],
        EXTENT*c(0.4, 0.45, 0.45, 0.4),
        border=2, col="#ff000033", lty=1)
    }
    par(op)
  })
  output$table_rem <- renderTable({
    m()$removal
  }, rownames = TRUE, colnames = TRUE, digits = 0)
  output$table_vis <- renderTable({
    m()$visits
  }, rownames = TRUE, colnames = TRUE, digits = 0)
  output$plot_est <- renderPlot({
    op <- par(mar=c(0,0,0,0))
    plot(m())
    par(op)
    print(str(e()))
  })
}

shinyApp(ui = ui, server = server)
