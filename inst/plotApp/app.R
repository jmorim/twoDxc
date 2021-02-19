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

ui = fluidPage(
  theme = shinytheme('darkly'),

  # User input panel -----------------------------------------------------------
  sidebarPanel(
    tabsetPanel(
      # file settings tab ------------------------------------------------------
      tabPanel(
        title = 'File',
        # file input -----------------------------------------------------------
        fileInput(
          inputId = 'filename',
          label = 'Data file(s)',
          multiple = T,
          accept = 'mzml'
        ),
        # file select ----------------------------------------------------------
        checkboxGroupInput(
          inputId = 'file.select',
          label = 'Select files to plot',
          choiceNames = NULL,
          choiceValues = NULL,
          selected = NULL
        ),
        # mod time -------------------------------------------------------------
        numericInput(
          inputId = 'mod.time',
          label = 'Modulation time (s)',
          value = 60
        ),
        # delay time -----------------------------------------------------------
        numericInput(
          inputId = 'delay.time',
          label = 'Delay time (s)',
          value = 30
        )
      ),
      # Plot params tab --------------------------------------------------------
      tabPanel(
        title = 'Plot Params',
        # plot type ------------------------------------------------------------
        selectInput(
          inputId = 'plot.type',
          label = 'Plot type',
          choices = c('2D', '3D')
        ),
        # rt1 range ------------------------------------------------------------
        sliderInput(
          inputId = 'rt1.range',
          label = '1D RT range (min)',
          value = c(0, 60),
          min = 0, max = 60,
          round = -1
        ),
        # rt2 range ------------------------------------------------------------
        sliderInput(
          inputId = 'rt2.range',
          label = '2D RT range (s)',
          value = c(0, 60),
          min = 0, max = 60
        ),
        # m/z range ------------------------------------------------------------
        sliderInput(
          inputId = 'mz.range',
          label = 'm/z range',
          value = c(100, 1700),
          min = 50, max = 1700,
          round = -4
        ),
        # eic ion and tolerance ------------------------------------------------
        checkboxInput(
          inputId = 'eic',
          label = 'EIC'
        ),
        fluidRow(
          # eic ion-------------------------------------------------------------
          column(7,
                 numericInput(
                   inputId = 'eic.ion',
                   label = 'Ion',
                   step = 0.0001,
                   value = NULL
                 ),
                 style = 'padding:0px;'
                 ),
          # eic tolerance ------------------------------------------------------
          column(5,
                 numericInput(
                   inputId = 'mz.tol',
                   label = 'ppm',
                   value = 25
                 ),
                 style = 'padding:0px;'
                 )
        ),
        # intensity range ------------------------------------------------------
        sliderInput(
          inputId = 'int.range',
          label = 'Intensity range',
          value = c(0, 0),
          min = 0, max = 0
        ),
        # color scale ----------------------------------------------------------
        selectInput(
          inputId = 'color',
          label = 'Color scale',
          choices = rownames(RColorBrewer::brewer.pal.info),
          selected = 'Spectral'
        ),
        # reverse scale --------------------------------------------------------
        checkboxInput(
          inputId = 'color.reverse',
          label = 'Reverse scale',
          value = T
        )
      )
    ),
    # update button ------------------------------------------------------------
    actionButton(
      inputId = 'update',
      label = 'Update'
    ),
    width = 3
  ),
  # Plot output panel ----------------------------------------------------------
  mainPanel(
    # make output take up 80% of window height
    tags$head(tags$style('.shiny-plot-output{height:80vh !important;}')),
    # tabs for output
    tabsetPanel(
      # xcms plot tab ----------------------------------------------------------
      tabPanel(
        title = 'xcms plot',
        plotOutput(
          outputId = 'xcms.plot'
          ) %>% withSpinner()
      ),
      # 2d plot tab ------------------------------------------------------------
      tabPanel(
        title = '2D plot',
        plotlyOutput(
          outputId = 'plotly.plot'
        ) %>% withSpinner()
      ),
      # test table tab ---------------------------------------------------------
      tabPanel(
        title = 'MSnObject table',
        tableOutput(
          outputId = 'test.table'
        )
      ),
      # calculators ------------------------------------------------------------
      tabPanel(
        title = 'MS calculators',
        fluidRow(
          # eic ion-------------------------------------------------------------
          column(3,
                 numericInput(
                   inputId = 'eic.ion.calc',
                   label = 'Ion',
                   step = 0.0001,
                   value = NULL
                 ),
                 style = 'padding:0px;'
                 ),
          # eic tolerance ------------------------------------------------------
          column(2,
                 numericInput(
                   inputId = 'mz.tol.calc',
                   label = 'ppm',
                   value = 25
                 ),
                 style = 'padding:0px;'
                 ),
          column(4,
                 selectInput(
                   inputId = 'adduct',
                   label = 'Adduct',
                   choices = c('M+H', 'M+Na', 'M+K', '2M+H'),
                   selected = 'M+H'
                 ))
        ),
        tableOutput(
          outputId = 'calc.table'
        ),
        tableOutput(
          outputId = 'adduct.table'
        )
      )
    )
  )
)
#===============================================================================
server = function(input, output, session){
  options(shiny.maxRequestSize = 1024^3,
          shiny.reactlog = T)
# Handle inputs ----------------------------------------------------------------
  # Get filenames from filename input ------------------------------------------
  file.data= reactive({
    req(input$filename)
    input$filename
  })
  # Read ms.data ---------------------------------------------------------------
  ms.data = reactive({
    readMSData(file.data()$datapath, mode = 'onDisk')
  })

  # Observer for updating file selection list ----------------------------------
  observe({
    req(file.data())

    updateCheckboxGroupInput(
      session,
      inputId = 'file.select',
      choiceNames = file.data()$name,
      choiceValues = file.data()$datapath,
      selected = file.data()$datapath
    )
  })
  # Indices of selected files --------------------------------------------------
  file.indices = reactive({
    req(file.data(), input$file.select)
    which(file.data()$datapath %in% input$file.select)
  }) %>% debounce(1100)

  # Observer for updating slider limits ----------------------------------------
  # change to observer instead of observeEvent?
  observeEvent(file.indices(),{
    # rt1.range ----------------------------------------------------------------
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
      max = rt.max
    )
    # rt2.range ----------------------------------------------------------------
    rt2.max = input$mod.time
    updateSliderInput(
      session,
      inputId = 'rt2.range',
      value = c(0, rt2.max),
      max = rt2.max
    )
    # m/z range ----------------------------------------------------------------
    min.mzs = c()
    max.mzs = c()
    for(i in file.indices()){
      min.mz = floor(min(
        filterFile(ms.data(), file = i)@featureData@data$lowMZ))
      max.mz = ceiling(max(
        filterFile(ms.data(), file = i)@featureData@data$highMZ))
      min.mzs = c(min.mzs, min.mz)
      max.mzs = c(max.mzs, max.mz)
    }
    mz.min = min(min.mzs)
    mz.max = max(max.mzs)

    updateSliderInput(
      session,
      inputId = 'mz.range',
      value = c(mz.min, mz.max),
      min = mz.min, max = mz.max
    )
    # intensity range ----------------------------------------------------------
    max.ints = c()
    for(i in file.indices()){
      max.int = ceiling(max(tic(
        filterFile(ms.data(), file = i))))
      max.ints = c(max.ints, max.int)
    }
    int.max = max(max.ints)

    updateSliderInput(
      session,
      inputId = 'int.range',
      value = c(0, int.max),
      max = int.max
    )
  })
  # Reactive for plot parameters -----------------------------------------------
  plot.params = eventReactive(input$update, {
    req(ms.data())

    # If eic checked, use eic ion, else ion is null
    isolate({
      if(input$eic == T){
        plot.eic.ion = input$eic.ion
      }else{
        plot.eic.ion = NULL
      }
    })

    # Change plot.type selection from 2D to 2Di for plotly
    plot.type = input$plot.type
    if(plot.type == '2D'){
      plot.type = '2Di'
    }

    # list params, returned to plot.params
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
      int.max = input$int.range[2] + 1,
      plot.type = plot.type
    )
  })

# Handle outputs ---------------------------------------------------------------
  # xcms plot ------------------------------------------------------------------
  output$xcms.plot = renderPlot({
    req(plot.params(), file.indices())
    # depend on update button

    plot.parameters = plot.params()

    n.files = length(file.indices())

    mz.range = isolate({
      if(input$eic == T && !is.null(plot.parameters$ion)){
        calc.mz.Window(plot.parameters$ion, plot.parameters$ppm.tol)
      }else{
        c(plot.params()$mz.min,
          plot.params()$mz.max)
      }
    })

    chrom = suppressMessages(
      chromatogram(
        object = filterFile(
          ms.data(), file = file.indices()
        ),
        rt = c(plot.params()$rt.min,
               plot.params()$rt.max),
        mz = mz.range,
        aggregationFun = 'sum'
      )
    )

    plot.colors = brewer.pal(n.files, plot.parameters$color.scale)
    names(plot.colors) = file.data()[file.indices(), 1]

    plot(
      chrom, col = plot.colors[file.indices()])
    legend('top', legend = file.data()[file.indices(), 1],
           col = plot.colors, pch = 16)
  })
  # 2d plot output -------------------------------------------------------------
  output$plotly.plot = renderPlotly({
    req(ms.data(), plot.params())

#    plot.parameters = plot.params()

    do.call(plot2D, plot.params())

#    plot2D(
#      object = plot.parameters$object,
#      file = 1,
#      mod.time = plot.parameters$mod.time,
#      delay.time = plot.parameters$delay.time,
#      ion = plot.parameters$ion,
#      ppm.tol = plot.parameters$ppm.tol,
#      color.scale = plot.parameters$color.scale,
#      reverse.scale = plot.parameters$reverse.scale,
#      rt.min = plot.parameters$rt.min,
#      rt.max = plot.parameters$rt.max,
#      rt2.min = plot.parameters$rt2.min,
#      rt2.max = plot.parameters$rt2.max,
#      mz.min = plot.parameters$mz.min,
#      mz.max = plot.parameters$mz.max,
#      int.min = plot.parameters$int.min,
#      int.max = plot.parameters$int.max,
#      plot.type = plot.parameters$plot.type
#      )
  })
  # test table output -----------------------------------------------------------
  output$test.table = renderTable({
    req(ms.data())

    ms.data = filterFile(ms.data(), file.indices())
    data = data.frame(rt = rtime(ms.data),
                      tic = tic(ms.data))
    head(data)
  })

  # calc table output ----------------------------------------------------------
  output$calc.table = renderTable({
    req(input$eic.ion.calc, input$mz.tol.calc)

    mz.window = calc.mz.Window(input$eic.ion.calc, input$mz.tol.calc)

    mz.diff = mz.window[2] - input$eic.ion.calc

    mz.table = c(input$eic.ion.calc, mz.window, mz.diff)
    mz.table = t(mz.table)
    colnames(mz.table) = c('ion', 'low m/z', 'high m/z', 'abs diff')
    mz.table
  },
  digits = 4)

  output$adduct.table = renderTable({
    req(input$eic.ion.calc, input$adduct)

#    mass.table = data.frame(H = 1.007276,
#                            Na = 22.989218,
#                            K = 38.963158)
    mass.table = c(1.007276, 22.989218, 38.963158)
    names(mass.table) = c('H', 'Na', 'K')

    mw = case_when(
      input$adduct == 'M+H' ~ input$eic.ion.calc - mass.table['H'],
      input$adduct == 'M+Na' ~ input$eic.ion.calc - mass.table['Na'],
      input$adduct == 'M+K' ~ input$eic.ion.calc - mass.table['K'],
      input$adduct == '2M+H' ~ (input$eic.ion.calc - mass.table['H']) / 2
    )

    adducts = data.frame(Adduct = c('M', 'M+H', 'M+Na', 'M+K', '2M+H'),
                         'Accurate Mass' = c(mw,
                                             mw + mass.table,
                                             mw * 2 + mass.table['H']))

    adducts = adducts %>%
      mutate('low m/z' = Accurate.Mass - (Accurate.Mass * 25 / 10^6),
             'high m/z' = Accurate.Mass + (Accurate.Mass * 25 / 10^6))
  },
  digits = 4)
}

#thematic_shiny()
shinyApp(ui = ui, server = server)
