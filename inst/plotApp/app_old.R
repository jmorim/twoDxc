library(shiny)
library(plotly)
library(MSnbase)
library(twoDxc)
library(dplyr)
library(ggplot2)
library(shinycssloaders)
library(shinythemes)
library(thematic)
library(RColorBrewer)

##################
# Inputs and outputs
# Inputs:
#   data file - check
#   mod time
#   delay time
#   2d or 3d - check
#   rt1 range - check
#   rt2 range - check
#   mz range - check
#   eic ion and tolerance - check
#   int range - check
#   color scale - check
# Outputs:
#   plot
#   x-axis, y-axis, z-axis
#
##################

ui <- fluidPage(
  theme = shinytheme('darkly'),
  sidebarPanel(
    tabsetPanel(
      # *Input() functions,
      tabPanel(
        title = "File",
        # file
        fileInput(
          inputId = "filename",
          label = "Data file(s)",
          multiple = T,
          accept = ".mzml"
        ),
        uiOutput(
          outputId = 'file.select.rui'
        ),
        # mod time
        numericInput(
          inputId = "mod.time",
          label = "Modulation time (s)",
          value = 60
        ),
        # delay time
        numericInput(
          inputId = "delay.time",
          label = "Delay time (s)",
          value = 30
        )
      ),
      # input plot params ------------------------------------------------------
      tabPanel(
        title = "Plot Params",
        # 2d or 3d
        selectInput(
          inputId = "dimensions",
          label = "Plot dimensions",
          #          choices = c("2D", "2Di", "3D")
          choices = c("2Di", "3D")
        ),
        # rt1 range
        sliderInput(
          inputId = "rt1.range",
          label = "1D RT range (min)",
          value = c(0, 60),
          min = 0, max = 60,
          round = -1
        ),
        # rt2 range
        sliderInput(
          inputId = "rt2.range",
          label = "2D RT range (s)",
          value = c(0, 60),
          min = 0, max = 60
        ),
        # m/z range, isolate
        sliderInput(
          inputId = "mz.range",
          label = "mz range",
          value = c(100, 1700),
          min = 50, max = 1700,
          round = -4
        ),
        # eic ion and tolerance
        checkboxInput(
          inputId = "eic",
          label = "EIC"
        ),
        fluidRow(
          column(7,
            numericInput(
              inputId = "eic.ion",
              label = "Ion",
              step = 0.0001,
              value = NULL
            ),
            style = "padding:0px;"
          ),
          column(5,
            numericInput(
              inputId = "mz.tol",
              label = "ppm",
              value = 25
            ),
            style = "padding:0px;"
          )
        ),

        # int range
        sliderInput(
          inputId = "int.range",
          label = "intensity range",
          value = c(0, 0),
          min = 0, max = 0
        ),
        #        uiOutput('int.range.control'),

        # color scale
        selectInput(
          inputId = "color",
          label = "Color scale",
          choices = rownames(RColorBrewer::brewer.pal.info),
          selected = "Spectral"
        ),
        # reverse scale
        checkboxInput(
          inputId = "color.reverse",
          label = "Reverse scale",
          value = T
        )
      )
    ),
    # For updating
    actionButton(
      inputId = "update",
      label = "Update"
    ),
    width = 3
  ),
  mainPanel(
    tags$head(tags$style(".shiny-plot-output{height:80vh !important;}")),
    tableOutput(outputId = "mz.window"),
    # plot
    tabsetPanel(
      tabPanel(
        title = "2D plot",
        plotlyOutput(outputId = "plot") %>% withSpinner(color = "#0dc5c1")
      ),
      tabPanel(
        title = "xcms plot",
        #    plotlyOutput(outputId = 'plotly_plot'),
        plotOutput(outputId = "xcms_plot") %>% withSpinner(color = "#0dc5c1"),
      )
    ),
    width = 9
  )
)


server <- function(input, output, session) {
  options(shiny.maxRequestSize = 1024^3,
          shiny.reactlog = T)

  ms.data = reactive({
    req(input$filename)

    readMSData(input$filename$datapath, mode= 'onDisk')
  })

  file.names = reactive({
    req(ms.data())
    input$filename
  })

  # file.select-----------------------------------------------------------------
  output$file.select.rui = renderUI({
    req(ms.data())
    checkboxGroupInput(
      inputId = 'file.select',
      label = 'Select files to plot',
      choiceNames = file.names()$name,
      choiceValues = file.names()$datapath,
      selected = file.names()$datapath
    )
  })

  file.indices = reactive({
    req(file.names(), input$file.select)
    which(file.names()[, 4] %in% input$file.select)
  })


  # update inputs based on files selected --------------------------------------
  observeEvent(ms.data(),{
    req(file.indices())
    # rt1 range ----------------------------------------------------------------
#    rt.max = round(max(rtime(ms.data())) / 60, 3) + 0.001
    browser()
    max.rts = c()
    for(i in file.indices()){
      max.rt = round(max(
        rtime(filterFile(ms.data(), file = i))
        ) / 60, 3) + 0.001
      max.rts = c(max.rts, max.rt)
    }
    rt.max = max(max.rts)

    updateSliderInput(
      session,
      inputId = 'rt1.range',
      value = c(0, rt.max),
      min = 0, max = rt.max
    )
    # rt2 range ----------------------------------------------------------------
    rt2.max = input$mod.time
    updateSliderInput(
      session,
      inputId = 'rt2.range',
      value = c(0, rt2.max),
      min = 0, max = rt2.max
    )
    # m/z range ----------------------------------------------------------------
    mz.min <- floor(min(ms.data()@featureData@data$lowMZ))
    mz.max <- ceiling(max(ms.data()@featureData@data$highMZ))
    updateSliderInput(session, "mz.range",
      value = c(mz.min, mz.max),
      min = mz.min, max = mz.max
    )
    # intensity range ----------------------------------------------------------
    int.max <- max(ms.data()@featureData@data$totIonCurrent) + 1
    updateSliderInput(session, "int.range",
      value = c(0, int.max),
      min = 0, max = int.max
    )
  })

  plot.params <- eventReactive(input$update, {
#    if (!is.null(input$data.files)) {
#      ms.data <- readMSData(input$data.files$datapath, mode = "onDisk")

      req(input$filename)
      if (input$eic == T) {
        plot.eic.ion <- input$eic.ion
      } else {
        plot.eic.ion <- NULL
      }
      #    }
      list(
        object = ms.data(),
        file = 1,
        mod.time = input$mod.time,
        delay.time = input$delay.time,
        ion = plot.eic.ion,
        ppm.tol = input$mz.tol,
        color.scale = input$color,
        reverse.scale = input$color.reverse,
        rt.min = input$rt1.range[1] * 60,
        rt.max = input$rt1.range[2] * 60,
        rt2.min = input$rt2.range[1],
        rt2.max = input$rt2.range[2],
        mz.min = input$mz.range[1],
        mz.max = input$mz.range[2],
        int.min = input$int.range[1],
        int.max = input$int.range[2],
        plot.type = input$dimensions
      )
#    }
  })

# 2D/3D plotly output ----------------------------------------------------------

  output$plot <- renderPlotly({
    req(input$filename, plot.params())#, input$update)#, nrow(input$filename) == 1)
    do.call(plot2D, plot.params())
  })

# xcms plot output -------------------------------------------------------------
  output$xcms_plot = renderPlot({
    req(ms.data(), file.indices())

    input$update

#    n.files = nrow(input$filename)
    n.files = length(file.indices())

    mz.range = isolate({
      if(input$eic == T && !is.null(input$eic.ion)){
        calc.mz.Window(input$eic.ion, input$mz.tol)
      }else{
        input$mz.range
      }
    })

    chrom = suppressMessages(
      chromatogram(
      object = filterFile(
        ms.data(), file = file.indices()
        ),
      rt = isolate(input$rt1.range * 60),
      mz = mz.range,
      aggregationFun = 'sum'
      )
    )

    plot.colors = brewer.pal(n.files, input$color)
    names(plot.colors) = file.names()[file.indices(), 1]

    plot(
      chrom,
      col = plot.colors)

#    legend(x = 1, y = 95,
#           legend = input$filename[, 1],
#           col = brewer.pal(n.files, input$color))
    legend('top', legend = file.names()[file.indices(), 1],
           col = plot.colors, pch = 15)
  })

# m/z window table output ------------------------------------------------------

  output$mz.window <- renderTable(
    {
      if (!is.null(input$eic.ion)) {
        mz.range <- calc.mz.Window(input$eic.ion, input$mz.tol)
        mz.range <- t(mz.range)
        colnames(mz.range) <- c("low mz", "high mz")
        mz.range
      }
      else {
        mz.range <- data.frame("low mz" = NA, "high mz" = NA)
        colnames(mz.range) <- c("low mz", "high mz")
        mz.range
      }
    },
    digits = 4
  )
}

#old junk that used to be in plot output
######
#    if (!is.null(input$data.files) && input$update > 0) {
#      if (input$dimensions == "2D") {
#        do.call(plot2D, plot.params())
#      } else {
#        do.call(plot2D, plot.params())
#      }
#    }
#
#    output$plotly_plot = renderPlotly({
#      if (!is.null(input$data.files) && input$update > 0) {
#        if (input$dimensions == "3D" || input$dimensions == "2Di"){
#          do.call(plot2D, plot.params())
#        }else{
#          do.call(plot2D, plot.params())
#        }
#      }

thematic_shiny()
shinyApp(ui = ui, server = server)
