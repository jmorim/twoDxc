extractIonChromatogram = function(object, file, mz.min, mz.max){
      suppressWarnings(
        eic.data <- as.data.frame(filterMz(filterFile(object, file),
                                          c(mz.min, mz.max)))
      )
      stopifnot("No ions found in this range"= nrow(eic.data) > 0)
      # need to change column i to 'intensity' to match TIC method
      colnames(eic.data)[4] = 'intensity'

      # smush intensities if multiple ions are within tolerancecs per rt
      if(any(duplicated(eic.data$rt)) == T){
        eic.data = eic.data %>%
          group_by(rt) %>%
          summarise(intensity = sum(intensity), .groups = 'drop')
      }

      # Make data.frame from rt values
#      eic.df = as.data.frame(cbind(rt, rt.2d))

      # Get intensity vector for filtered ions an join to rt data, filling in rt
      # where there's 0 int
      eic.data = eic.data %>%
#        full_join(eic.data, by = 'rt') %>%
#        dplyr::select(rt, rt.2d, intensity)
        dplyr::select(rt, intensity)
      return(eic.data)
}

addScanRate <- function(x){
  stopifnot(inherits(x, 'MSnExp'))
  rts.by.file = split(rtime(x), fromFile(x))
  nScans = lapply(rts.by.file, length)
  x@experimentData@other$scanRate = lapply(seq_along(nScans), function(x){
    nScans[[x]] / (max(rts.by.file[[x]]) - min(rts.by.file[[x]]))
  })
  x
}

rasterize = function(object,
                     file = 1,
                     mod.time,
                     delay.time = 0,
                     ion = NULL,
                     mz.tol = c('ppm', 'abs'),
                     ppm.tol = 20,
                     abs.tol = 0.5,
                     mz.min = object@featureData$lowMZ,
                     mz.max = object@featureData$highMZ){

  # Handle args
  if(class(object) != 'MSnExp' & class(object) != 'OnDiskMSnExp'){
    stop('Error: not an MSnExp object')
  }
  if(!is.null(ion)){
    if(mz.tol != 'ppm' & mz.tol != 'abs'){
      stop('Error: Specify ppm or abs for m/z tolerance')
    }
  }

  # Get rt vector from data file, includes all retention times
  rt <- rtime(filterFile(object, file))
  # Convert these to 2D rts
  rt.2d = sapply(rt, convert.2drt, mod.time = mod.time,
                 delay.time = delay.time)

  # If not an EIC
  if(is.null(ion)){
    # but if mz range specified
    if(mz.min != object@featureData@data$lowMZ ||
       mz.max != object@featureData@data$highMZ){
      eic.data = extractIonChromatogram(object, file, mz.min, mz.max)
      # Join intensity to rt data
      plot.2d.df = as.data.frame(cbind(rt, rt.2d))
      plot.2d.df = plot.2d.df %>%
        full_join(eic.data, by = 'rt')
    # TIC case
    }else{
      # Get total ion current
      intensity = object@featureData@data %>%
        filter(fileIdx == file) %>%
        pull(totIonCurrent)
      # Join intensity to rt data
      plot.2d.df = as.data.frame(cbind(rt, rt.2d, intensity))
    }
    # EIC case
  }else{
    # if tol is ppm
    if(mz.tol == 'ppm'){
      mz.range = calc.mz.Window(ion, ppm.tol)
    # else tol is abs
    }else{
      mz.range = ion + c(-abs.tol, abs.tol)
    }
    eic.data = extractIonChromatogram(object,
                                       file,
                                       mz.range[1], mz.range[2])
      # Join intensity to rt data
    plot.2d.df = as.data.frame(cbind(rt, rt.2d))
    plot.2d.df = plot.2d.df %>%
      full_join(eic.data, by = 'rt')
  }

  # Group 1D retention times based on mod time
  rt.adj = NULL
  # raster number starts at delay time
  raster.counter = delay.time
  for(i in 1:length(rt)){
    # Set rt.adj to counter
    rt.adj[i] = raster.counter
    # If this rt is greater than the last one (e.g. 59.77 vs 0.02) bump raster
    if(rt.2d[i] > rt.2d[i + 1] & i != length(rt)){
      raster.counter = raster.counter + mod.time
    }
  }

  # Group 2D retention times so there aren't empty spots
  object = addScanRate(object)
  scan.rate = object@experimentData@other$scanRate[[file]]
  # Bin using scan rate
  # Need to test other scan rates to see if this holds
  rt.2d.adj = round(rt.2d * round(scan.rate)) / round(scan.rate)
#  rt.2d.adj = round(rt.2d * 26.5) / 26.5
#  data.points = length(rt)
#  scans.per.mod = ceiling(data.points / (max(rt) / mod.time))
#  rt.2d.adj = seq(1 / scan.rate, mod.time, length.out = scans.per.mod)
#  remainder = data.points - (length(unique(rt.adj)) * length(rt.2d.adj))
#
#  browser()

  # Add intensity, replace NAs with 0
  plot.2d.df$intensity[is.na(plot.2d.df$intensity)] = 0

  # Bind everything together
  plot.2d.df = cbind(plot.2d.df, rt.adj, rt.2d.adj)

  return(plot.2d.df)
}

#' Plot a 2D Chromatogram
#'
#' This function takes an MSnExp object and plots a 2D chromatogram using the
#' specified modulation time and optional ion extraction. The function relies on
#' ggplot2 for plotting.
#'
#' @param object An MSnExp object
#' @param file Which file (numeric, not filename) to pull data from. Default is
#'               1
#' @param mod.time The modulation time of the 2D run
#' @param delay.time The time it takes for modulation to start.
#'          Default is 0.
#' @param ion Select an ion to generate an extracted ion chromatogram (EIC). If
#'          missing, the total ion chromatogram (TIC) is returned.
#' @param mz.tol Whether EIC m/z range is calculated using ppm (relative) or
#'          absolute tolerance. Default is ppm
#' @param ppm.tol PPM tolerance, used to calculate m/z range for EIC. Default is
#'          20
#' @param abs.tol Absolute tolerance, used to calculate m/z range for EIC.
#'          Default is 0.5 m/z
#' @param log.scale Whether intensity coloring should be log scaled. Default is
#'          FALSE
#' @param rt.min Minimum 1D RT to plot. Default is 0
#' @param rt.max Maximum 1D RT to plot. Default is the max found in the data
#' @param rt2.min Minimum 2D RT to plot. Default is 0
#' @param rt2.max Maximum 2D RT to plot. Default is the modulation time
#' @param mz.min Minimum m/z to plot. Default is lowest found in the data
#' @param mz.max Maximum m/z to plot. Default is highest found in the data
#' @param int.min Minimum intensity to plot. Default is 0
#' @param int.max Maximum intensity to plot. Default is highest found in the TIC
#' @param color.scale Color scale from distiller to use. Default is 'Spectral.'
#' @param reverse.scale Reverse the color scale. Default is T.
#' @param save.output Whether to save the output as a file. Default is FALSE
#' @param filename The name for saving the file if desired. Default is
#'          '2dplot.png' plus the file selected plus the ion selected.
#' @param filetype Filetype of the form ".pdf". Default is ".pdf"
#' @param filepath Where to save the output file. Default is '.'
#' @param print.output Whether to show the plot output. Default is TRUE
#' @param mz.digits The number of digits after the decimal to show for m/z's
#'          in the plot title. Default is 4
#' @param plot.type Choose from 2D ggplot, interactive 2D plotly plot ('2Di'),
#'          or an interactive 3D plotly plot ('3D'). Default is 2D
#' @return A 2-dimensional plot
#' @importFrom reshape2 dcast
#' @import plotly
#' @export
setGeneric('plot2D', function(object,
                              file = 1,
                              mod.time,
                              delay.time = 0,
                              ion = NULL,
                              mz.tol = c('ppm', 'abs'),
                              ppm.tol = 20,
                              abs.tol = 0.5,
                              log.scale = F,
                              rt.min = 0,
                              rt.max = max(rtime(object)),
                              rt2.min = 0,
                              rt2.max = mod.time,
                              mz.min = min(object@featureData$lowMZ),
                              mz.max = max(object@featureData$highMZ),
                              int.min = 0,
                              int.max = max(tic(object), na.rm = T)+1,
                              color.scale = 'Spectral',
                              reverse.scale = T,
                              save.output = F,
                              filename = '2dplot',
                              filetype = '.pdf',
                              filepath = '.',
                              print.output = T,
                              mz.digits = 4,
                              plot.type = c('2D', '2Di', '3D'))
           standardGeneric('plot2D'))

setMethod('plot2D', 'MSnExp', function(object,
                                       file = 1,
                                       mod.time,
                                       delay.time = 0,
                                       ion = NULL,
                                       mz.tol = c('ppm', 'abs'),
                                       ppm.tol = 20,
                                       abs.tol = 0.5,
                                       log.scale = F,
                                       rt.min = 0,
                                       rt.max = max(rtime(object)) + delay.time,
                                       rt2.min = 0,
                                       rt2.max = mod.time,
                                       mz.min = min(object@featureData$lowMZ),
                                       mz.max = max(object@featureData$highMZ),
                                       int.min = 0,
                                       int.max = ceiling(
                                         max(
                                           filterFile(
                                             object, file = file
                                           )@featureData@data$totIonCurrent
                                           )
                                         )+ 1,
#                                           tic(object), na.rm = T)) + 1,
                                       color.scale = 'Spectral',
                                       reverse.scale = T,
                                       save.output = F,
                                       filename = '2dplot',
                                       filetype = '.pdf',
                                       filepath = '.',
                                       print.output = T,
                                       mz.digits = 4,
                                       plot.type = c('2Di', '2D', '3D')){

  mz.tol = mz.tol[1]
  plot.type = plot.type[1]

  if(!is.null(ion)){
    mz.range = calc.mz.Window(ion, ppm.tol)
    mz.min = mz.range[1]
    mz.max = mz.range[2]
  }

  plot.data = rasterize(
    object,
    file = file,
    mod.time = mod.time,
    delay.time = delay.time,
    ion = ion,
    mz.tol = mz.tol,
    ppm.tol = ppm.tol,
    abs.tol = abs.tol,
    mz.min = mz.min,
    mz.max = mz.max
  )

  # Change intensity to log scale if specified
  if(log.scale == T){
    plot.data$intensity[is.na(plot.data$intensity)] = 1
    plot.data = plot.data %>%
      mutate(intensity = log(intensity))
  }

  # Trim data according to specified limits
  plot.data = plot.data %>%
    # rt lim
    filter(between(rt.adj, rt.min, rt.max)) %>%
    # rt2 lim
    filter(between(rt.2d.adj, rt2.min, rt2.max)) %>%
    # intensity
    filter(between(intensity, int.min, int.max))

  # Check desired return format. If wide, needs to be formatted so that
  # columns represent rt1, rows represent rt2, and values represent intensity
  # Useful for plotly plots
  if(plot.type == '2Di' || plot.type == '3D'){
    plot.data = reshape2::dcast(plot.data[,3:5],
                       rt.2d.adj ~ rt.adj,
                       value.var = 'intensity',
                       fun.aggregate = mean)
    rownames(plot.data) = plot.data$rt.2d.adj
    plot.data = plot.data[, -1]
  }
  # Turn scale direction into number, T == 1, F == 0
  scale.direction = 1 - (2 * reverse.scale)

  if(plot.type== '2D'){
    plot.2d = ggplot(plot.data, aes(x = rt.adj, y = rt.2d.adj)) +
      geom_raster(aes(fill = intensity), interpolate = T, na.rm = T) +
      scale_fill_distiller(palette = color.scale,
                           direction = scale.direction) +
      xlab('1D Retention Time (s)') +
      ylab('2D Retention Time (s)') +
      ggtitle(
        if(is.null(ion)){
          paste0('TIC: File ', file)
        }else{
          paste0('EIC: Ions ', round(mz.min, digits = mz.digits), ' -- ',
                 round(mz.max, digits = mz.digits), ', File: ', file)
        }) +
      theme_classic()
    if(save.output == T){
      if(is.null(ion) & filename == '2dplot'){
        ggsave(paste0(filename, filetype), plot = plot.2d, path = filepath)
      }else{
        filename = paste0(filename, '_file_', file, '_', ion, filetype)
        ggsave(paste0(filename), plot = plot.2d, path = filepath)
      }
    }
  }else if(plot.type== '2Di'){
    plot.2d = plot_ly(z = as.matrix(plot.data),
                      type = 'heatmap',
                      colors = color.scale,
                      reversescale = reverse.scale,
                      zsmooth = 'best') %>%
      layout(xaxis = list(title = '1D Retention Time (min)'),
             yaxis = list(title = '2D Retention Time (s)'))
  }else{
    plot.2d = plot_ly(z = ~as.matrix(plot.data),
                      colors = color.scale,
                      reversescale = reverse.scale) %>%
      add_surface()
  }
  plot.2d
})
