library(shiny)

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

ui = fluidPage(
  # *Input() functions,
  # file
  fileInput(inputId = 'data.files',
            label = 'Data file(s)',
            multiple = T,
            accept = '.mzml'),
  # mod time
  numericInput(inputId = 'mod.time',
               label = 'Modulation time (s)',
               value = 60),
  # delay time
  numericInput(inputId = 'delay.time',
               label = 'Delay time (s)',
               value = 30),
  # 2d or 3d
  selectInput(inputId = '2d.3d',
              label = 'Plot dimensions',
              choices = c('2D', '3D')),
  # rt1 range
  sliderInput(inputId = 'rt1.range',
              label = '1D RT range',
              value = c(1, 60),
              min = 1, max = 60),
  # rt2 range
  sliderInput(inputId = 'rt2.range',
              label = '2D RT range',
              value = c(1, 60),
              min = 1, max = 60),
  # m/z range
  sliderInput(inputId = 'mz.range',
              label = 'mz range',
              value = c(50, 1000),
              min = 50, max = 1000),
  # eic ion and tolerance
  checkboxInput(inputId = 'eic',
                label = 'EIC'),
  numericInput(inputId = 'eic.ion',
               label = 'EIC ion',
               value = NULL),
  numericInput(inputId = 'mz.tol',
               label = 'm/z tolerance (ppm)',
               value = 25),

  # int range
  sliderInput(inputId = 'int.range',
              label = 'intensity range',
              value = c(0, 10^9),
              min = 50, max = 10^9),

  # color scale
  selectInput(inputId = 'color',
              label = 'Color scale',
              choices = rownames(RColorBrewer::brewer.pal.info),
              selected = 'Spectral'),
  # reverse scale
  checkboxInput(inputId = 'color.reverse',
                label = 'Reverse scale'),

  # *Output() functions,
  # mz window
  tableOutput(outputId = 'mz.window'),

  # plot
  plotOutput(outputId = 'plot')
)

server = function(input, output) {
  options(shiny.maxRequestSize = 1024*1024^2)
#  ms.data = readMSData(input$data.files, mode = 'onDisk')
#  output$plot = renderPlot({
#    library(MSnbase)
#    ms.data = readMSData(input$data.files, mode = 'onDisk')
#
#    plot2D()
#  })
  output$mz.window = renderTable({
    if(!is.null(input$eic.ion)){
      mz.range = calc.mz.Window(input$eic.ion, input$mz.tol)
      mz.range = t(mz.range)
      colnames(mz.range) = c('low mz', 'high mz')
      mz.range
    }
    else{
      mz.range = data.frame('low mz' = NA, 'high mz' = NA)
      colnames(mz.range) = c('low mz', 'high mz')
      mz.range
    }
  },
#  colnames = T,
  digits = 4)
}

calc.mz.Window <- function(mz, ppm){
  lower.mz <- mz - (mz * ppm / 10^6)
  upper.mz <- mz + (mz * ppm / 10^6)
  return(c(lower.mz, upper.mz))
}

shinyApp(ui = ui, server = server)
