# TODO
# - Add 2D RT constraints for filtering matching psgs
#   - DONE Write helper function to calculate 2D RT from mod.time and dead.time
#   - Fix rtmin.2d and rtmax.2d because right now max can be less than min
# - INPROG Make calc.mz.Window optionally take absolute mz tolerances, not just ppm
# - Add argument for minimum number of qualifier ions to group pseudospectra
#    together
# - Add regression testing argument to filter out bad ions after 2D grouping
# - DONE Add function that plots 2D chromatogram heatmap
#     - Limit digits for EIC in plot title, default to 4 past 0
# - Remove empty elements from string columns in pspec.2d, e.g. adducts
# - DONE Count the 2D peak groups
# - Allow subsetting of 2D plots to minimize processing time

setClass('xsAnnotate2D',
         representation(pspec2D = 'list',
                        mod.time = 'numeric',
                        dead.time = 'numeric'),
         contains=c('xsAnnotate'))

# Make group2D function
#' Group 2D Peaks
#'
#' This function takes an xsAnnotate object produced by CAMERA and groups the
#' pseudospectrum groups according to their m/z's, 1D retention time, and 2D
#' retention time.
#' The default chromatographic/MS tolerances are 5s 2D RT tolerance and 20 ppm
#' m/z tolerance. Your system's mileage may vary.
#'
#' @param object An xsAnnotate object returned by CAMERA that needs 2D peak
#'          grouping
#' @param mod.time The modulation time for the comprehensive 2D run
#' @param dead.time The time before the injection reaches the detector, defaults
#'          to 0
#' @param rt.tol The grouping tolerance in the RT dimension
#' @param ppm.tol The grouping tolerance in the m/z dimension
#' @param parallelized If true, runs grouping in parallel. Requires the future
#'          package
#' @return An xsAnnotate2D object with the added pspec2D slot
#' @export
setGeneric('group2D', function(object, mod.time, dead.time = 0, rt.tol = 5,
                                ppm.tol = 20, parallelized = F)
  standardGeneric('group2D'))
setMethod('group2D', 'xsAnnotate', function(object, mod.time, dead.time = 0,
                                             rt.tol = 5, ppm.tol = 20,
                                             parallelized = F){
  if (class(object) != 'xsAnnotate'){
    stop('Error: not an xsAnnotate object')
  }


  # Function for matching mzs with apply loop
  # Not sure if it's bad to put a function in a function
  matchMzs <- function(grouped.psg.mz, all.pspecs, ppm.tol, rt.tol){
      # Calculate the m/z range with ppm.tol
      mz.range = calc.mz.Window(grouped.psg.mz, ppm.tol)
      # Get the 1DRT

      # Calculate the 2drt range with rt.tol
#      rt.range =
      # Get matching m/z's within tolerance in the other pspectra
      # Add in RT tolerance
      matching.mzs = all.pspecs %>%
        filter(mz > mz.range[1] & mz < mz.range[2])
      # Initialize matrix for condensed ion data
      condensed.ion = NULL
      # Populate matrix with summary functions to squish m/z data together
      condensed.ion <- matching.mzs %>%
        # mz col based on mean mz
        summarise(mz = mean(mz)) %>%
        bind_cols(matching.mzs %>%
                    # take corresponding mins and maxes
                    summarise(mzmin = min(mzmin)),
                  matching.mzs %>%
                    summarise(mzmax = max(mzmax)),
                  # Take RT from most intense peak modulation (max in 1D direction)
                  matching.mzs %>%
                    summarise(rt = matching.mzs[which(
                      getInt(matching.mzs) ==
                        max(getInt(matching.mzs), na.rm = T),
                      arr.ind = T)[1], 'rt']),
                  matching.mzs %>%
                    summarise(rtmin = min(rtmin)),
                  matching.mzs %>%
                    summarise(rtmax = max(rtmax)))
      # sum all peaks from npeaks to (but not including) isotopes,
      # including intensities
      # Need to adapt this for single sample case
      if("npeaks" %in% colnames(matching.mzs)){
        condensed.ion <- condensed.ion %>%
          bind_cols(matching.mzs %>%
                      select(npeaks:isotopes) %>%
                      select(-isotopes) %>%
                      summarise_all(sum, na.rm = T))
      }else{
        condensed.ion <- condensed.ion %>%
          bind_cols(matching.mzs %>%
                      summarise_at('into', sum))
      }
      # not a perfect solution for summarizing string data but can
      # improve later
      condensed.ion <- condensed.ion %>%
        bind_cols(matching.mzs %>%
                    summarise(isotopes = list(isotopes),
                    adducts = list(adduct),
                    psgs = list(psg))) %>%
      # calculate 2d rts
        mutate(rt.2d = convert.2drt(rt, mod.time, dead.time),
             rtmin.2d = convert.2drt(rtmin, mod.time, dead.time),
             rtmax.2d = convert.2drt(rtmax, mod.time, dead.time)) %>%
        # Remove redundant ion data
        distinct(mz, rt, .keep_all = T) %>%
        # Reorder data
        select(mz:rtmax,rt.2d:rtmax.2d, everything())
    #  mz.counter <<- mz.counter + 1
      return(condensed.ion)
    }

  # Function for matching psgs with matched m/zs and/or rts
  matchPsgs <- function(pseudospec, all.pspecs, ppm.tol, parallelized = F){
    # Get pspectrum data for one pspec
    pspecgroup = all.pspecs %>%
      filter(psg == pseudospec)
    # Take average retention time
    #rt = mean(pspecgroup[, 'rt'])
    # Determine quantitative ion by finding most common max ion in each sample
    max.ion.vector = apply(getInt(pspecgroup), 2, which.max)
    # Get the mode max ion's index, use to get the most intense ion
    max.ion = pspecgroup[getMode(max.ion.vector), 'mz']
    # Calculate m/z tolerance window using ppm.tol
    mz.window = calc.mz.Window(max.ion, ppm.tol)
    # Find psgs with ions in this m/z window
    matching.psgs <- all.pspecs %>%
      filter(mz > mz.window[1] & mz < mz.window[2]) %>%
      # pull out psg number
      pull(psg)
    # Put all matching psgs into one 2D peak group
    grouped.psgs <- all.pspecs %>%
      filter(psg %in% matching.psgs)
    condensed.psg <- NULL
    # List all m/z's
    mzs <- grouped.psgs[, 'mz']
    # Find matching m/z's and add intensities together using matchMzs, above
    if(parallelized == T){
      condensed.psg <- do.call(rbind, future_lapply(
        mzs, matchMzs, all.pspecs = all.pspecs, ppm.tol = ppm.tol,
        rt.tol = rt.tol))
    }else{
    condensed.psg <- do.call(rbind, lapply(
      mzs, matchMzs, all.pspecs = all.pspecs, ppm.tol = ppm.tol,
      rt.tol = rt.tol))
    }
    # Remove redundant ions from the psg
    condensed.psg <- condensed.psg %>%
      distinct(mz, rt, .keep_all = T) %>%
      mutate('psg.2d' = psg.counter)
    # Print output message
    cat(paste0('Adding 2D pspec ', psg.counter, '.\n'))
    # Add to counter
    psg.counter <<- psg.counter + 1
    return(condensed.psg)
  }

  # Convert xsAnnotate to xsAnnotate2D
  new.object <- as(object, 'xsAnnotate2D')
  new.object@mod.time = mod.time
  new.object@dead.time = dead.time

  # Get all pseudospectra (pspecs) stored in the object
  all.pspectra <- CAMERA::getpspectra(new.object,
                                      grp = 1:length(object@pspectra))
  # List all the pspec groups (psgs)
  psgs <- unique(all.pspectra[, 'psg'])

  ### Debugging length, comment when done
  psgs <- psgs[1:10]

  # Initialize 2D pspec list and global variables
  pspec2D <- list()
  # Don't really need an m/z counter at the moment
 # mz.counter <<- 1
  psg.counter <- 1
#  mod.time <- mod.time
#  dead.time <- dead.time

  # Check if future.apply is loaded if running in parallel is specified.
  # For some reason !requireNamespace doesn't ever return TRUE
#  if(parallelized == T && !requireNamespace('future.apply', quietly = T)){
  if(parallelized == T && !'future.apply' %in% loadedNamespaces()){
    warning("The package 'future.apply' must be installed to run group2D in
            parallel. Defaulting to serial processing.\n")
    parallelized = F
  }
  # For each psg, use the function matchPsgs, see below
  if(parallelized == T){
    pspec2D <- do.call(rbind, future_lapply(
      psgs, matchPsgs, all.pspecs = all.pspectra, ppm.tol = 20,
      parallelized = T))
    # Might change to using furrr and purrr
#    pspec2D <- do.call(rbind, future_map(
#      psgs, matchPsgs, all.pspecs = all.pspectra, ppm.tol = 20,
#      parallelized = T))
  }else{
    pspec2D <- do.call(rbind, lapply(
      psgs, matchPsgs, all.pspecs = all.pspectra, ppm.tol = 20))
    # Might change to using furrr and purrr
#    pspec2D <- do.call(rbind, map(
#      psgs, matchPsgs, all.pspecs = all.pspectra, ppm.tol = 20))
  }

  # Remove redundant pspecs (ones with same m/z and rt)
  # This is only problematic with co-eluting isomers, in which case there's
  # nothing you can do anyway.
  cat('Condensing...\n')
  pspec2D.distinct <- pspec2D %>%
    distinct(mz, rt, .keep_all = T)

  # The psg.2d number assignments are going to be staggered because duplicates
  # are removed, so to change e.g. 1, 1, 1, 3, 3, 5, to 1, 1, 1, 2, 2, 3:
  old.pspec2D <- pspec2D.distinct
  j <- 1
  for(i in 1:length(pspec2D.distinct$psg.2d)){
    if(i > 1){
      if(pspec2D.distinct$psg.2d[i] != old.pspec2D$psg.2d[i - 1]){
        j <- j + 1
      }
      pspec2D.distinct$psg.2d[i] = j
    }
  }

  # Calculate 2D retention time and reorder
#  pspec2D.distinct <- pspec2D.distinct %>%
#    mutate(rt.2d = convert.2drt(rt, mod.time, dead.time),
#           rtmin.2d = convert.2drt(rtmin, mod.time, dead.time),
#           rtmax.2d = convert.2drt(rtmax, mod.time, dead.time)) %>%
#    select(mz:rtmax,rt.2d:rtmax.2d, everything())

  # Set the object's pspec2D slot
  new.object@pspec2D <- pspec2D.distinct

  # Print output message showing how many pspecs were grouped
  cat('Grouped', length(unique(new.object@pspec2D$psg.2d)), '2D pseudospectra\n')
  # Remove global vars
#  rm('psg.counter', 'mod.time', 'dead.time', pos = 1L)
  # Return object
  return(new.object)
})

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
#' @param dead.time The time it takes for the injection to reach the detector.
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
#' @param save.output Whether to save the output as a file. Default is FALSE
#' @param filename The name for saving the file if desired. Default is
#'          '2dplot.png' plus the file selected plus the ion selected.
#' @param filepath Where to save the output file. Default is '.'
#' @param print.output Whether to show the plot output. Default is TRUE
#' @param mz.digits The number of digits after the decimal to show for m/z's
#'          in the plot title. Default is 4
#' @return A 2-dimensional plot
#' @export
setGeneric('plot2D', function(object, file = 1, mod.time, dead.time = 0, ion,
                              mz.tol = c('ppm', 'abs'), ppm.tol = 20,
                              abs.tol = 0.5, log.scale = F,
                              save.output = F, filename = '2dplot.png',
                              filepath = '.', print.output = T, mz.digits = 4)
  standardGeneric('plot2D'))

setMethod('plot2D', 'MSnExp', function(object, file = 1, mod.time,
                                       dead.time = 0,
                                       ion, mz.tol = c('ppm', 'abs'),
                                       ppm.tol = 20, abs.tol = 0.5,
                                       log.scale = F,
                                       save.output = F, filename = '2dplot.png',
                                       filepath = '.',
                                       print.output = T, mz.digits = 4){
  ###
  # Is there a better way to select first element in list as default arg?
  mz.tol <- mz.tol[1]
  ###
  if(class(object) != 'MSnExp' & class(object) != 'OnDiskMSnExp'){
    stop('Error: not an MSnExp object')
  }
  if(missing(ion)){
    if(mz.tol != 'ppm' & mz.tol != 'abs'){
      stop('Error: Specify ppm or abs for m/z tolerance')
    }
  }

  # Get rt vector from data file, includes all retention times
  rt <- filterFile(object, file)@featureData@data$retentionTime
  # Convert these to 2D rts
  rt.2d <- sapply(rt, convert.2drt, mod.time = mod.time, dead.time = dead.time)
  if(missing(ion)){
    # Get total ion current or extracted ion current if ion specified
    intensity <- object@featureData@data %>%
      filter(fileIdx == file) %>%
      pull(totIonCurrent)
    # Join intensity to rts data
    plot.2d.df <- as.data.frame(cbind(rt, rt.2d, intensity))
  }else{
    if(mz.tol == 'ppm'){
      mz.range <- calc.mz.Window(ion, ppm.tol)
    }else{
      mz.range <- ion + c(-abs.tol, abs.tol)
    }
    # filterMz data to get EIC data of structure [file, rt, mz, i]
    # Suppress warnings because each call to trimMz_Spectrum from filterMz
    # is very helpful and tells me when a scan has an empty spectrum.
    suppressWarnings(
      eic.data <- as.data.frame(filterMz(filterFile(object, file), mz.range))
    )
    # Need to change column name 'i' to 'intensity' to match TIC method
    colnames(eic.data)[4] = 'intensity'

    # Need to smush intensities in cases where there are >1 m/z per rt
    if(any(duplicated(eic.data$rt)) == T){
      eic.data <- eic.data %>%
        group_by(rt) %>%
        # hopefully this works because it seems to return a goofy table
        summarise(intensity = sum(intensity))
    }
    # Make data.frame from rt values like in TIC version
    plot.2d.df <- as.data.frame(cbind(rt, rt.2d))
    # Get just the intensity vector for the filtered ions
    # Join intensity to rts data, filling in rts where there's 0 int
    plot.2d.df <- plot.2d.df %>%
      full_join(eic.data, by = 'rt') %>%
      select(rt, rt.2d, intensity)
  }

  # Group 1D retention times based on modulation time
  rt.adj <- NULL
  raster.counter <- 0
  for(i in 1:length(rt)){
    if(rt.2d[i] > rt.2d[i + 1] & i != length(rt)){
      raster.counter <- raster.counter + mod.time
    }
    rt.adj[i] = raster.counter
  }

  # Group 2D retention times so there aren't any empty spots
  rt.2d.adj = round(rt.2d)

  # Use log scale if specified (useful for EICs) and remove NAs
  # Need to use 1 if log scale so -Inf doesn't show up
  if(log.scale == T){
    plot.2d.df$intensity[is.na(plot.2d.df$intensity)] <- 1
    plot.2d.df <- plot.2d.df %>%
      mutate(intensity = log(intensity))
  }else{
    plot.2d.df$intensity[is.na(plot.2d.df$intensity)] <- 0
  }

  # plot function
    plot.2d <- ggplot(plot.2d.df, aes(x = rt.adj, y = rt.2d.adj)) +
      geom_raster(aes(fill = intensity), interpolate = T) +
      scale_fill_distiller(palette = 'Spectral', direction = -1) +
      xlab('1D Retention Time (s)') +
      ylab('2D Retention Time (s)') +
      ggtitle(
        if(missing(ion)){
          paste0('TIC: File ', file)
        }else{
          paste0('EIC: Ions ', round(mz.range[1], digits = mz.digits), ' - ',
                 round(mz.range[2], digits = mz.digits), ', File: ', file)
      }) +
      theme_classic()
  if(save.output == T){
    if(missing(ion) & filename == '2dplot.png'){
      ggsave(filename, plot = plot.2d, path = filepath)
    }else{
      filename <- paste0(unlist(strsplit(filename, '.png')), '_file_', file,
                         '_', ion, '.png')
      ggsave(filename, plot = plot.2d, path = filepath)
    }
  }
  if(print.output == T){
    print(plot.2d)
  }
})

# Function to get intensities from pseudospectrum
# Need to modify to take intensities from single sample datasets
getInt <- function(x){
  if("npeaks" %in% colnames(x)){
    int.mat <- x %>%
      select(matches('X\\d'))
  }else{
    int.mat <- x %>%
      select(into)
  }
  return(int.mat)
}

# Function to get mode from a vector
# Used to pick the quantitative ion in a pseudospectrum
getMode <- function(x){
  u <- unique(x)
  mode <- u[which.max(tabulate(match(x, u)))]
  mode <- unlist(mode)
  return(mode)
}

# Get quantitative (main) ion by intensity
                  # Take RT from most intense peak modulation (max in 1D direction)
#                  matching.mzs %>%
#                    summarise(rt = matching.mzs[which(
#                      getInt(matching.mzs) ==
#                        max(getInt(matching.mzs), na.rm = T),
#                      arr.ind = T)[1], 'rt'])
getQuantIon <- function(x){
  q.ion <- x[which(getInt(x) == max(getInt(x), na.rm = T), arr.ind = T)[1],
            'mz']
  return(q.ion)
}

#' @export
# Function to calculate mz range given mz and ppm tolerance
# Returns range in a list
calc.mz.Window <- function(mz, ppm){
  lower.mz <- mz - (mz * ppm / 10^6)
  upper.mz <- mz + (mz * ppm / 10^6)
  return(c(lower.mz, upper.mz))
}

#' @export
# Function to calculate 2D RT from the 1D RT
convert.2drt <- function(rt, mod.time, dead.time = 0) {
  rt.adj <- rt - dead.time
  rt.2d <- rt.adj - (mod.time * floor(rt.adj / mod.time))
  return(rt.2d)
}
