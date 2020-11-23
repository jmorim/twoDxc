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
                        delay.time = 'numeric'),
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
#' @param delay.time The time before the injection reaches the detector, defaults
#'          to 0
#' @param rt.tol The grouping tolerance in the RT dimension
#' @param ppm.tol The grouping tolerance in the m/z dimension
#' @param parallelized If true, runs grouping in parallel. Requires the future
#'          package
#' @param tail.fix If true, corrects for erroneous grouping when prevalent
#'   compounds tail in the first dimension
#'
#' @return An xsAnnotate2D object with the added pspec2D slot
#' @export
setGeneric('group2D', function(object, mod.time, delay.time = 0, rt.tol = 180,
                               rt2.tol = 5, ppm.tol = 20, parallelized = F,
                               tail.fix = F, verbose = F)
  standardGeneric('group2D'))
setMethod('group2D', 'xsAnnotate', function(object, mod.time, delay.time = 0,
                                             rt.tol = 180, rt2.tol = 5,
                                            ppm.tol = 20, parallelized = F,
                                            tail.fix = F, verbose = F){
  if (class(object) != 'xsAnnotate'){
    stop('Error: not an xsAnnotate object')
  }

  # Convert xsAnnotate to xsAnnotate2D
  new.object <- as(object, 'xsAnnotate2D')
  new.object@mod.time <- mod.time
  new.object@delay.time <- delay.time

  # Get all pseudospectra (pspecs) stored in the object
  all.pspectra <- CAMERA::getpspectra(new.object,
                                      grp = 1:length(object@pspectra))
  # List all the pspec groups (psgs)
  psgs <- unique(all.pspectra[, 'psg'])

  ### Debugging length, comment when done
#  psgs <- psgs[1:10]

  # Initialize 2D pspec list and global variables
  pspec2D <- list()
  # Don't really need an m/z counter at the moment
 # mz.counter <<- 1
  psg.counter <<- 1L
#  psg.counter <- 1L
#  count.env <- new.env()
#  local(psg.counter <- 1L, env = count.env)
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
  ########
  # Parallel procecssing for psgs not supported yet because I don't know
  # how to set counters for parallel processes
#  if(parallelized == T){
#    pspec2D <- do.call(rbind, future_lapply(
#      psgs, matchPsgs, all.pspecs = all.pspectra, ppm.tol = ppm.tol,
#      parallelized = T))
    # Might change to using furrr and purrr
##    pspec2D <- do.call(rbind, future_map(
##      psgs, matchPsgs, all.pspecs = all.pspectra, ppm.tol = 20,
##      parallelized = T))
#  }else{
    pspec2D <- do.call(rbind, lapply(
      psgs, matchPsgs, all.pspecs = all.pspectra, ppm.tol = ppm.tol,
      parallelized = parallelized, verbose = verbose,
      mod.time = mod.time, delay.time = delay.time, rt.tol = rt.tol,
      rt2.tol = rt2.tol))
    # Might change to using furrr and purrr
#    pspec2D <- do.call(rbind, map(
#      psgs, matchPsgs, all.pspecs = all.pspectra, ppm.tol = 20))
#  }

  # Remove redundant pspecs (ones with same m/z and rt)
  # This is only problematic with co-eluting isomers, in which case there's
  # nothing you can do anyway.
  cat('Added ', nrow(pspec2D), ' spectra\n')
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

#  # Print output message showing how many pspecs were grouped
#  cat('Grouped', length(unique(new.object@pspec2D$psg.2d)), '2D pseudospectra\n')
  # Remove global vars, find a better way to do this
  rm(psg.counter, pos = 1L)

  if(tail.fix == T){
    new.pspecs <- data.frame()
    for(i in 1:nrow(new.object@pspec2D)){
      mz = new.object@pspec2D[i, 'mz']
      rt = new.object@pspec2D[i, 'rt']
      rt.2d = new.object@pspec2D[i, 'rt.2d']
      x.vec = c(mz, rt, rt.2d)
      # Look for rows within mz, rt, and rt.2 limits
      dup.psgs <-new.object@pspec2D %>%
        filter(near(mz, x.vec[1], tol = x.vec[1] * ppm.tol / 10^6),
               near(rt, x.vec[2], tol = rt.tol),
               near(rt.2d, x.vec[3], tol = rt2.tol))
      # Pick the integration with max signal
      psg.int <- getInt(dup.psgs)
      total.psg <- psg.int %>%
        replace(is.na(.), 0) %>%
        mutate(total = rowSums(.))
      max.psg <- which.max(total.psg$total)
      # Add max psg integration to new df
      new.pspecs <- rbind(new.pspecs, dup.psgs[max.psg,])
    }
    new.pspecs <- new.pspecs %>%
      distinct(mz, rt, rt.2d, .keep_all = T)
    new.object@pspec2D <- new.pspecs
    # Print output message showing how many pspecs were grouped
    cat('Grouped', length(unique(new.object@pspec2D$psg.2d)), '2D pseudospectra\n')
    return(new.object)
  }else{
  # Print output message showing how many pspecs were grouped
  cat('Grouped', length(unique(new.object@pspec2D$psg.2d)), '2D pseudospectra\n')
  # Return object
  return(new.object)
  }
})

# Function for matching mzs with apply loop
# Not sure if it's bad to put a function in a function
matchMzs <- function(mzRt, grouped.psg.rt, all.pspecs,
                     ppm.tol, rt.tol, rt2.tol,
                     mod.time = mod.time, delay.time = delay.time){
    # Get m/z and RT data
  range.195 = calc.mz.Window(195.0880, 25)
    mz <- mzRt[1]
#    if(mz > range.195[1] & mz < range.195[2]){
#      browser()
#    }
    rt <- mzRt[2]
    # Calculate the m/z range with ppm.tol
    mz.range <- calc.mz.Window(mz, ppm.tol)
    # Calculate the RT ranges with rt.tol and rt2.tol
    rt.range <- rt + c(-rt.tol, rt.tol)
    rt2.range <- convert.2drt(rt, mod.time = mod.time,
                              delay.time = delay.time) + c(-rt2.tol, rt2.tol)

    # Get matching m/z's within tolerance in the other pspectra
    # Add in RT tolerance
    matching.mzs = all.pspecs %>%
      filter(mz > mz.range[1] & mz < mz.range[2]) %>%
      filter(rt > rt.range[1] & rt < rt.range[2]) %>%
      mutate(rt.2d = convert.2drt(rt, mod.time, delay.time),
             rtmin.2d = convert.2drt(rtmin, mod.time, delay.time),
             rtmax.2d = convert.2drt(rtmax, mod.time, delay.time)) %>%
      filter(rt.2d > rt2.range[1] & rt.2d < rt2.range[2])
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
                # Take RT from most intense peak modulation
                # (max in 1D direction)
                matching.mzs %>%
                  summarise(rt = matching.mzs[which(
                    getInt(matching.mzs) ==
                      max(getInt(matching.mzs), na.rm = T),
                    arr.ind = T)[1], 'rt']),
                matching.mzs %>%
                  summarise(rtmin = min(rtmin)),
                matching.mzs %>%
                  summarise(rtmax = max(rtmax)),
                matching.mzs %>%
                  summarise(rt.2d = matching.mzs[which(
                    getInt(matching.mzs) == max(getInt(matching.mzs),
                                                na.rm = T),
                    arr.ind = T)[1], 'rt.2d']),
                matching.mzs %>%
                  summarise(rtmin.2d = min(rtmin.2d)),
                matching.mzs %>%
                  summarise(rtmax.2d = max(rtmax.2d))
      )
    # sum all peaks from npeaks to (but not including) isotopes,
    # including intensities
    # Need to adapt this for single sample case
    if("npeaks" %in% colnames(matching.mzs)){
      condensed.ion <- condensed.ion %>%
        bind_cols(matching.mzs %>%
                    select(npeaks:isotopes) %>%
                    select(-isotopes) %>%
                    select(-ms_level) %>%
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
      # Remove redundant ion data
      distinct(mz, rt, .keep_all = T) %>%
      # Reorder data
      select(mz:rtmax,rt.2d:rtmax.2d, everything())
  #  mz.counter <<- mz.counter + 1
    return(condensed.ion)
}

# Function for matching psgs with matched m/zs and/or rts
matchPsgs <- function(pseudospec, all.pspecs, ppm.tol, rt.tol, rt2.tol,
                      parallelized = F, verbose = F,
                      mod.time = mod.time, delay.time = delay.time){
  #Sys.sleep(1)
#  if(psg.counter == 21){
#    browser()
#  }
#   browser()
  # Get pspectrum data for one pspec
  pspecgroup <- all.pspecs %>%
    filter(psg == pseudospec)

  # Determine quantitative ion by finding most common max ion in each sample
  max.ion.vector <- apply(getInt(pspecgroup), 2, which.max)
  # Get the mode max ion's index, use to get the most intense ion
  max.ion <- pspecgroup[getMode(max.ion.vector), 'mz']
  # Calculate m/z tolerance window using ppm.tol
  mz.window <- calc.mz.Window(max.ion, ppm.tol)

  # Get quant ion's RT
  max.ion.rt <- pspecgroup[getMode(max.ion.vector), 'rt']
  # Calculate RT ranges with rt.tol and rt2.tol
  rt.range <- max.ion.rt + c(-rt.tol, rt.tol)
  rt2.range <- convert.2drt(max.ion.rt, mod.time = mod.time,
                            delay.time = delay.time) + c(-rt2.tol, rt2.tol)

  # Find psgs with ions in this m/z and RT window
  matching.psgs <- all.pspecs %>%
    filter(mz > mz.window[1] & mz < mz.window[2]) %>%
    filter(rt > rt.range[1] & rt < rt.range[2]) %>%
    mutate(rt.2d = convert.2drt(rt, mod.time, delay.time),
           rtmin.2d = convert.2drt(rtmin, mod.time, delay.time),
           rtmax.2d = convert.2drt(rtmax, mod.time, delay.time)) %>%
    filter(rt.2d > rt2.range[1] & rt.2d < rt2.range[2]) %>%
    # pull out psg number
    pull(psg)
  # Put all matching psgs into one 2D peak group
  grouped.psgs <- all.pspecs %>%
    filter(psg %in% matching.psgs)
  condensed.psg <- NULL
  # List all m/z's
  mzs <- grouped.psgs[, 'mz']
  # Figure out a way to condense this list, no need to run matchMzs for
  # 195.0878 and 195.0877
  rts <- grouped.psgs[, 'rt']
  mzRts <- grouped.psgs %>%
    select(mz, rt)
  # Find matching m/z's and add intensities together using matchMzs, above
  if(parallelized == T){
    condensed.psg <- do.call(rbind, future_apply(
      mzRts, MARGIN = 1, matchMzs, all.pspecs = all.pspecs, ppm.tol = ppm.tol,
      rt.tol = rt.tol, rt2.tol = rt2.tol,
      mod.time = mod.time, delay.time = delay.time))
  }else{
  condensed.psg <- do.call(rbind, apply(
    mzRts, MARGIN = 1, matchMzs, all.pspecs = all.pspecs, ppm.tol = ppm.tol,
    rt.tol = rt.tol, rt2.tol = rt2.tol,
    mod.time = mod.time, delay.time = delay.time))
#    condensed.psg <- do.call(rbind, pmap(
#      mzRts, matchMzs, all.pspecs = all.pspecs, ppm.tol = ppm.tol,
#      rt.tol = rt.tol, rt2.tol = rt2.tol,
#      mod.time = mod.time, dead.time = dead.time
#    ))
  }
  # Remove redundant ions from the psg
  #browser()
#  psg.counter.here <- local(psg.counter, count.env)
#  psg.counter.here <- psg.counter
  condensed.psg <- condensed.psg %>%
    distinct(mz, rt, .keep_all = T) %>%
    mutate('psg.2d' = psg.counter)
  # Print output message
  if(verbose == T){
    cat(paste0('Adding 2D pspec ', psg.counter, '.\n'))
  }
  # Add to counter
  ## Fix this, check if there's a counter already in the df and add to it
  ## should get around the parallel problem.
  ## can't access pspec.2D from here, it won't be formed till apply is done
  psg.counter <<- psg.counter + 1L
#  local(psg.counter <- psg.counter.here + 1L, count.env)
  return(condensed.psg)
}

# Subtract ion function
#' Subtract ion
#'
#' This function subtracts an ion from the MSnbase object by removing its
#' intensity from the TIC.
#'
#' @param object An xsAnnotate object returned by CAMERA that needs 2D peak
#'          grouping
#' @param ion The ion to subtract
#' @param tol The tolerance window for the ion in ppm
#'
#' @return An xsAnnotate2D object with the ion intensity subtracted
#' @export
subtractIon = function(object, ion, tol){
  mz.range = calc.mz.Window(ion, tol)
  eic.df = extractIonChromatogram(object, file = 1, mz.range[1], mz.range[2])
  #  object.tic =
  tic.data = data.frame(rt = rtime(object), intensity = tic(object))
  subtracted.data = tic.data %>%
    full_join(eic.df, by = 'rt', suffix = c('.tic', '.eic')) %>%
    mutate(intensity = intensity.tic - intensity.eic) %>%
    pull(intensity)
  object@featureData@data$totIonCurrent = as.numeric(subtracted.data)
  return(object)
}

# Function to get intensities from pseudospectrum
# Need to modify to take intensities from single sample datasets
getInt <- function(x){
  if("npeaks" %in% colnames(x) | "npeaks" %in% names(x)){
    int.mat <- as.data.frame(x) %>%
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
  # Remove 0s
  u <- u[u!=0]
  mode <- u[which.max(tabulate(match(x, u)))]
  mode <- unlist(mode)
  return(mode)
}

# Get quantitative (main) ion by intensity
getQuantIon <- function(x){
  if(nrow(x) == 1){
    return(x[, 'mz'])
  }
  q.ion <- x[which(getInt(x) == max(getInt(x), na.rm = T), arr.ind = T)[1],
            'mz']
  return(q.ion)
}

# calc.mz.Window function
#' Function to calculate m/z window given an ion and ppm tolerance
#' @param mz m/z
#' @param ppm ppm tolerance for the m/z window
#'
#' @return A vector of length 2 with min m/z and max m/z
#' @export
calc.mz.Window <- function(mz, ppm){
  lower.mz <- mz - (mz * ppm / 10^6)
  upper.mz <- mz + (mz * ppm / 10^6)
  return(c(lower.mz, upper.mz))
}

# Function to calculate 2D RT from the 1D RT
convert.2drt <- function(rt, mod.time, delay.time = 0) {
  rt.adj <- rt - delay.time
  rt.2d <- rt.adj - (mod.time * floor(rt.adj / mod.time))
  return(rt.2d)
}
