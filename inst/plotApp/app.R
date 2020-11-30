library(shiny)
library(plotly)
library(MSnbase)
library(twoDxc)
library(dplyr)
library(ggplot2)

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
  sidebarPanel(
    tabsetPanel(
      # *Input() functions,
      tabPanel(
        title = "File",
        # file
        fileInput(
          inputId = "data.files",
          label = "Data file(s)",
          multiple = T,
          accept = ".mzml"
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
      tabPanel(
        title = "Plot",
        # 2d or 3d
        selectInput(
          inputId = "dimensions",
          label = "Plot dimensions",
          choices = c("2D", "2Di", "3D")
        ),
        # rt1 range
        sliderInput(
          inputId = "rt1.range",
          label = "1D RT range (min)",
          value = c(0, 60),
          min = 0, max = 72,
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
    tableOutput(outputId = "mz.window"),
    # plot
    plotOutput(outputId = "plot"),
#    plotlyOutput(outputId = 'plotly_plot'),
    width = 9
  )
)


server <- function(input, output, session) {
  options(shiny.maxRequestSize = 1024^3,
          shiny.reactlog = T)
  #  ms.data = readMSData(input$data.files, mode = 'onDisk')

  # update slider values based on input file
  #  output$int.range.control = renderUI({
  #    if(!is.null(input$data.files)){
  #      ms.data = readMSData(input$data.files$datapath, mode = 'onDisk')
  #      max.int = max(ms.data@featureData@data$totIonCurrent)
  #    sliderInput(inputId = 'int.range', label = 'Intensity range',
  #                min = 0, max = max.int, value = c(0, max.int))
  #    }
  #  })

  #  output$range.controls
  # change from observe to reactive?
  observe({
    if (!is.null(input$data.files)) {
      ms.data <- readMSData(input$data.files$datapath, mode = "onDisk")
    } else {
      ms.data <- readMSData(system.file(
        "tea_data", "prime.mzML",
        package = "twoDxc"
      ), mode = "onDisk")
    }
    # update controls based on data file
    # rt1
    rt.max <- round(max(ms.data@featureData@data$retentionTime) / 60, 3) + 0.001
    updateSliderInput(session, "rt1.range",
      value = c(0, rt.max),
      min = 0, max = rt.max
    )
    # rt2
    rt2.max <- input$mod.time
    updateSliderInput(session, "rt2.range",
      value = c(0, rt2.max),
      min = 0, max = rt2.max
    )
    # mz
    mz.min <- floor(min(ms.data@featureData@data$lowMZ))
    mz.max <- ceiling(max(ms.data@featureData@data$highMZ))
    updateSliderInput(session, "mz.range",
      value = c(mz.min, mz.max),
      min = mz.min, max = mz.max
    )
    # intensity
    int.max <- max(ms.data@featureData@data$totIonCurrent) + 1
    updateSliderInput(session, "int.range",
      value = c(0, int.max),
      min = 0, max = int.max
    )
  })

  plot.params <- eventReactive(input$update, {
    if (!is.null(input$data.files)) {
      ms.data <- readMSData(input$data.files$datapath, mode = "onDisk")
      if (input$eic == T) {
        plot.eic.ion <- input$eic.ion
      } else {
        plot.eic.ion <- NULL
      }
      #    }
      list(
        object = ms.data,
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
    }
  })

  output$plot <- renderPlot({
    if (!is.null(input$data.files) && input$update > 0) {
      if (input$dimensions == "2D") {
        do.call(plot2D, plot.params())
      } else {
        do.call(plot2D, plot.params())
      }
    }

#    output$plotly_plot = renderPlotly({
#      if (!is.null(input$data.files) && input$update > 0) {
#        if (input$dimensions == "3D" || input$dimensions == "2Di"){
#          do.call(plot2D, plot.params())
#        }else{
#          do.call(plot2D, plot.params())
#        }
#      }
    }
    )

    #    plot2D(plot.params())
    #    if(!is.null(input$data.files)){
    #      ms.data = readMSData(input$data.files$datapath, mode = 'onDisk')
    #      if(input$eic == T){
    #        plot.eic.ion = input$eic.ion
    #      }else{
    #        plot.eic.ion = NULL
    #      }
    #      plot2D(ms.data, file = 1, input$mod.time, input$delay.time,
    #             ion = plot.eic.ion, ppm.tol = input$mz.tol,
    #             color.scale = input$color, reverse.scale = input$color.reverse,
    #             rt.min = input$rt1.range[1] * 60, rt.max = input$rt1.range[2] * 60,
    #             mz.min = isolate({input$mz.range[1]}),
    #             mz.max = isolate({input$mz.range[2]}))
    #    }
    #
    #    plot2D()

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

# calc.mz.Window <- function(mz, ppm){
#  lower.mz <- mz - (mz * ppm / 10^6)
#  upper.mz <- mz + (mz * ppm / 10^6)
#  return(c(lower.mz, upper.mz))
# }

shinyApp(ui = ui, server = server)
