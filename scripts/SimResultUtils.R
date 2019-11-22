
#' Approximate harmonic mean of a continuous function by calculating the
#' mean across the function value at \code{steps} evenly-spaced points 
#' between \code{x1} and \code{x2}.
#' 
#' @param traj  The function to calculate the mean for.
#' @param x1,x2 The interval over which to calculate the mean.
#' @param steps The number of evenly-spaced points at which to 
#'              calculate the function value for approximating 
#'              the mean.
#'              
approxHarmonicMean <- function(traj, x1, x2, steps=100) {
  x <- seq(x1, x2, length.out=steps)
  return ( 1/mean( 1/traj(x) ) )
}


#' Approximate arithmetic mean of a continuous function by calculating the
#' mean across the function value at \code{steps} evenly-spaced points.
#' between \code{x1} and \code{x2}
#'
#' @inheritParams approxHarmonicMean
#' 
approxArithmeticMean <- function(traj, x1, x2, steps=100) {
  x <- seq(x1, x2, length.out=steps)
  return ( mean(traj(x)) )
}


#' The geometric mean of of the approximate harmonic and arithmetic means 
#' of a continuous function by calculating the mean across the function value
#' at \code{steps} evenly-spaced points between \code{x1} and \code{x2}.
#' 
#' \deqn{ \sqrt(\frac{a}{\Sum{f(x)^{-1}}} b \sum{f(x)} }
#' 
#' @inheritParams approxHarmonicMean
#' @param a,b The weighting for the harmonic and geometric means, respectively
#' 
approxGeomCombMean <- function(traj, x1, x2, a=1, b=1, steps=100) {
  x <- seq(x1, x2, length.out=steps)
  return(sqrt( (a/mean( 1/traj(x)) * b*mean(traj(x)))))
}


#' Approximate the range of a continuous function by calculating the range at 
#' \code{steps} evenly-spaced points between \code{x1} and \code{x2}.
#' 
#' @param traj  The function to calculate the range for.
#' @param x1,x2 The interval over which to calculate the range.
#' @param steps The number of evenly-spaced points at which to 
#'              calculate the function value for approximating 
#'              the range.
#'              
approxRange <- function(traj, x1, x2, steps=100) {
  x <- seq(x1, x2, length.out=steps)
  res <- range(traj(x))
  names(res) <- c("min", "max")
  return(res)
}


#' Return a step function of the mean, calculated using \code{meanFn}, 
#' between every two points in the timegrid \code{t}. Also return 
#' the range of the function within each interval.
getMeanTraj <- function(t, traj, meanFn, ...) {
  
  parTruth <- c()
  parRange <- c()
  for (j in 2:length(t)) {
    parTruth <- c(parTruth, meanFn(traj, t[j-1], t[j], ...))
    parRange <- rbind(parRange, approxRange(traj, t[j-1], t[j]))
  }
  
  return(list(truth=parTruth, range=parRange))
}


# Not used
piecewiseFn <- function(x, vals, breaks) {
  if (length(x) > 1) {
    return(sapply(x, piecewiseFn, vals, breaks))
  } else {
    return( vals[max(1, which.min(breaks < x) - 1)] )
  }
}

# Not used
getPiecewiseFn <- function(vals, breaks) {
  return (function(x) piecewiseFn(x, vals, breaks))
}




#' Get summary statistics for a single simulation replicate, by comparing the true trajectory simulated under to the 
#' median and HPD interval of the inference.
#'
#' The function will return a list with several summary statistics:
#' - bias:               Mean of the difference between the median posterior estimate and the truth, across each interval.
#' - rel.bias:           Bias divided by the truth.
#' - weighted.bias:      Mean of the difference between the median posterior estimate and the truth, across each interval and weighted by the length of each interval.
#' 
#' @param t The times of interval boundaries. Should start with 0 and be of dimension (n+1)
#' @param truthmean The mean (or harmonic mean, depending on the mean function used) of the true trajectory for each 
#'                  interval. Should be of dimension n
#' @param truthrange The min and max of the true trajectory within each interval. Should have column names (min,max) 
#'                   and have n rows.
#' @param hpdest The HPD and median estimates for the parameter in each interval. Should have column names (lower, med, upper)
#'               and have n rows.
#'
getRunStats <- function(t, truthmean, truthrange, hpdest) {
  
  if (!is.matrix(hpdest)) {
      hpdest <- as.matrix(t(hpdest))
  }
  
  # Mean bias, absolute deviation and mean credible interval width
  bias     <- mean(hpdest[, "med"] - truthmean)
  mad      <- mean(abs(hpdest[,"med"] - truthmean))
  mciw     <- mean(hpdest[, "upper"] - hpdest[, "lower"])
  
  # Mean relative absolute deviation and mean relative credible interval width
  rel.bias <- mean((hpdest[, "med"] - truthmean)/truthmean)
  rel.mad  <- mean(abs(hpdest[,"med"] - truthmean)/truthmean)
  rel.mciw <- mean((hpdest[, "upper"] - hpdest[, "lower"])/truthmean)
  
  # Weighted mean absolute deviation and weighted mean credible interval width
  weighted.bias <- weighted.mean(hpdest[, "med"] - truthmean, diff(t))
  weighted.mad  <- weighted.mean(abs(hpdest[,"med"] - truthmean), diff(t))
  weighted.mciw <- weighted.mean(hpdest[, "upper"] - hpdest[, "lower"], diff(t))
  
  # Weighted mean absolute deviation and weighted mean credible interval width
  rel.weighted.bias <- weighted.mean((hpdest[, "med"] - truthmean)/truthmean, diff(t))
  rel.weighted.mad  <- weighted.mean(abs(hpdest[,"med"] - truthmean)/truthmean, diff(t))
  rel.weighted.mciw <- weighted.mean((hpdest[, "upper"] - hpdest[, "lower"])/truthmean, diff(t))
  
  # Coverage of the true mean and coverage of the true range by the HPD
  meancov <- mean(truthmean >= hpdest[, "lower"] & truthmean <= hpdest[, "upper"])
  abscov  <- mean((truthrange[, "min"] >= hpdest[, "lower"] & truthrange[,"min"] <= hpdest[, "upper"]) & 
                    (truthrange[, "max"] >= hpdest[, "lower"] & truthrange[,"max"] <= hpdest[, "upper"]))
  
  # Weighted coverage of the true mean and weighted coverage of the true range by the HPD
  weighted.meancov <- weighted.mean(truthmean >= hpdest[, "lower"] & truthmean <= hpdest[, "upper"], diff(t))
  weighted.abscov  <- weighted.mean((truthrange[, "min"] >= hpdest[, "lower"] & truthrange[,"min"] <= hpdest[, "upper"]) & 
                                      (truthrange[, "max"] >= hpdest[, "lower"] & truthrange[,"max"] <= hpdest[, "upper"]), diff(t))
  
  return(data.frame(bias = bias,  rel.bias = rel.bias, weighted.bias = weighted.bias, rel.weighted.bias = rel.weighted.bias,
                    mad  = mad,   rel.mad  = rel.mad,  weighted.mad  = weighted.mad,  rel.weighted.mad  = rel.weighted.mad, 
                    mciw = mciw,  rel.mciw = rel.mciw, weighted.mciw = weighted.mciw, rel.weighted.mciw = rel.weighted.mciw, 
                    meancov = meancov, abscov = abscov, weighted.meancov = weighted.meancov, weighted.abscov = weighted.abscov))
}



#' Get the summary statistics for all replicates of a simulation scenario (trajectory)
#' 
#' @param simresult
#' @param popSizes
#' @param sampIntensities
#' @param popSizeTmies
#' @param sampIntensityTimes
#' @param cutPopSizes If TRUE population size summary statistics are only computed over the interval from the oldest to the most recent sample
#'                    (the interval from th tMRCA to the oldest sample is left out). If FALSE population size summary statistics are computed
#'                    from the tMRCA to the present. Sampling intensity summary statistics are always only computed over the interval from the 
#'                    oldest to the most recent sample, since the sampling intensity is not defined between the tMRCA and the oldest sample.
#' @param plotResults If TRUE plot four plots for each replicate in different directories, prefixed by outputPath:
#'                    - figures/:
#'                    - figures-log/:
#'                    - figures-short/:
#'                    - figures-short-log/:
#' @param outputPath
getStats <- function(simresult, popSizes, sampIntensities=NULL, popSizeTimes, sampIntensityTimes=NULL, cutPopSizes=TRUE, plotResults=TRUE, outputPath="") {
  
  if (outputPath != "" && plotResults == TRUE) {
      dir.create(paste0(outputPath,"figures/"), recursive=TRUE, showWarnings=FALSE)
      dir.create(paste0(outputPath,"figures-log/"), recursive=TRUE, showWarnings=FALSE)
      dir.create(paste0(outputPath,"figures-short/"), recursive=TRUE, showWarnings=FALSE)
      dir.create(paste0(outputPath,"figures-short-log/"), recursive=TRUE, showWarnings=FALSE)
      
      outnames <- sapply(strsplit(chanames(popSizes), "/"), function(x) x[length(x)])
      outnames <- sapply(strsplit(outnames, "\\."), function(x) paste(x[1:(length(x)-1)], collapse="."))
  }
  
  popSizeStats <- sampIntensityStats <- c()
  for (i in 1:length(popSizes)) {
      
      # PopSize HPDs and segment times
      popSizeHPD <- getHPDMedian(popSizes[[i]])
      t <- c(0, popSizeTimes[[i]])
      
      # Sampling intensity HPDs and epoch times
      maxSampleTime <- max(simresult[[i]]$samp_times)    
      if (!is.null(sampIntensities)) {
          sampIntensityHPD <- getHPDMedian(sampIntensities[[i]])
          s <- c(0, sampIntensityTimes[[i]])  
          
          # If the analysis had a density-defined sampling protocol
          if (length(sampIntensityTimes[[i]]) == 1) {
              samplingDim        <- length(simresult[[i]]$samp_intensity)
              sampIntensityHPD   <- matrix(rep(sampIntensityHPD, samplingDim), 
                                           nrow=samplingDim, byrow=TRUE, dimnames=list(1:samplingDim, names(sampIntensityHPD)))
              s <- c(simresult[[i]]$epoch_times[1:samplingDim], sampIntensityTimes[[i]])
          }
          
          # If the simulation had a density-defined sampling protocol
          if (simresult[[i]]$nrepochs == 1) {
              simresult[[i]]$samp_intensity <- rep(simresult[[i]]$samp_intensity, length(sampIntensityTimes[[i]]))
          }
          
          if (abs(maxSampleTime - s[length(s)]) > 1e-6) {
              stop(paste("Oldest sampling time mismatch between simulation and analysis:",maxSampleTime,"!=", s[length(s)]))
          }
      } else {
          sampIntensityHPD <- NULL
          sampIntensityTruth <- NULL
          s <- c(0, maxSampleTime)
      }
      
      
      # Disregard everything after the oldest sampling time
      if (cutPopSizes == TRUE && s[length(s)] < t[length(t)]) {
          j <- which.min(t < s[length(s)])
          t <- t[1:j]
          t[j] <- s[length(s)]
          popSizeHPD <- popSizeHPD[1:(j-1),]
      }
  
      # True trajectory
      traj    <- function(x) simresult[[i]]$traj(x + simresult[[i]]$timeshift)
      logtraj <- function(x) log( simresult[[i]]$traj(x + simresult[[i]]$timeshift) )
      
      # Harmonic mean of population size in each interval (continuous function)
      popSizeTruth <- c()
      popSizeRange <- c()
      for (j in 2:length(t)) {
        popSizeTruth <- c(popSizeTruth, approxHarmonicMean(traj, t[j-1], t[j]))
        popSizeRange <- rbind(popSizeRange, approxRange(traj, t[j-1], t[j]))
      }
      
      # Calculate statistics  
      popSizeStats       <- rbind(popSizeStats, getRunStats(t, popSizeTruth, popSizeRange, popSizeHPD))
      
      # True samping intensity in each interval (piecewise constant function => f(x) = mean(f(x)) = range(f(x)))
      if (!is.null(sampIntensities)) {
          sampIntensityTruth <- simresult[[i]]$samp_intensity
          sampIntensityRange <- cbind(sampIntensityTruth, sampIntensityTruth)
          colnames(sampIntensityRange) <- c("min","max")
          
          # Calculate statistics
          sampIntensityStats <- rbind(sampIntensityStats, getRunStats(s, sampIntensityTruth, sampIntensityRange, sampIntensityHPD))
      }
  
      
      # Plot replicate
      if (outputPath != "" && plotResults == TRUE) {
          plotSimResult(simresult[[i]]$phylo, 
                        popSizeHPD, t, traj, popSizeTruth, 
                        sampIntensityHPD, s, sampIntensityTruth, 
                        logPopSize=FALSE, short=FALSE, outputfile=paste0(outputPath, "figures/", outnames[i],".pdf"))
          
          plotSimResult(simresult[[i]]$phylo, 
                        popSizeHPD, t, traj, popSizeTruth, 
                        sampIntensityHPD, s, sampIntensityTruth, 
                        logPopSize=FALSE, short=TRUE, outputfile=paste0(outputPath, "figures-short/", outnames[i],".pdf")) 
          
          plotSimResult(simresult[[i]]$phylo, 
                        popSizeHPD, t, traj, popSizeTruth, 
                        sampIntensityHPD, s, sampIntensityTruth, 
                        logPopSize=TRUE, short=FALSE, outputfile=paste0(outputPath, "figures-log/", outnames[i],".pdf"))
          
          plotSimResult(simresult[[i]]$phylo, 
                        popSizeHPD, t, traj, popSizeTruth, 
                        sampIntensityHPD, s, sampIntensityTruth, 
                        logPopSize=TRUE, short=TRUE, outputfile=paste0(outputPath, "figures-short-log/", outnames[i],".pdf")) 
      }
  }
  
  if (outputPath != "") {
      if (cutPopSizes == TRUE) {
          write.csv(popSizeStats,       paste0(outputPath, "popSizeStatsSamplingPeriod.csv"), row.names=FALSE)
      } else {
          write.csv(popSizeStats,       paste0(outputPath, "popSizeStatsFull.csv"), row.names=FALSE)
      }
    
      if (!is.null(sampIntensities)) {
          write.csv(sampIntensityStats, paste0(outputPath, "sampIntensityStats.csv"), row.names=FALSE)
      }
  }
  
  if (!is.null(sampIntensities)) {
      return(list(popSizeStats = popSizeStats, sampIntensityStats = sampIntensityStats))
  } else {
      return(list(popSizeStats = popSizeStats))
  }
}







plotStats <- function(stats, ylab="", ylim=NULL, names=c(), plotGrid=TRUE, plotStrip=TRUE, las=1) {
  
  
  if (is.null(ylim)) {
      ymin <- min(c(0, sapply(stats, min)))
      ymin <- max(ymin, -2)
      
      ymax <- max(sapply(stats, max))
      ymax <- min(ymax, 2)
      
      ylim <- c(ymin, ymax)  
  }
  
  boxplot(stats, ylab = ylab, outline=ifelse(plotStrip, FALSE, TRUE), col=mPal(oxCols$blue1,0.5),
          ylim=ylim, names = names, las = las, lwd = 1, cex.lab=1.5, cex.axis=1.2, xaxs='i', yaxs='i', xaxt='n')
  
  if (plotGrid) {
      grid(col=mPal(oxCols$gray6), nx=0, ny=NULL)
      #abline(h=axTicks(2), lty=2, col=mPal(oxCols$gray3), lwd=0.5)
  }
  
  if (plotStrip) {
      stripchart(stats, vertical = TRUE, method = "jitter", jitter = 0.2, add = TRUE, pch = 20, col = mPal(oxCols$red2, 0.25), cex=1)
  }
  
}


plotStatsComparison <- function(stats, ylab="", ylim=NULL, names=c(), plotGrid=TRUE, plotStrip=FALSE, las=1, cols=mPal(oxCols$blue2)) {
  
  n <- length(stats)
  x <- 1 + (n+1)*(0:(length(stats[[1]])-1))
  
  cols <- rep(cols, n)
  
  boxplot(stats[[1]], ylab = ylab, outline=ifelse(plotStrip, FALSE, TRUE), col=cols[1],
          shrink=0.8, ylim=ylim, xlim=c(0,length(stats[[1]])*(n+1)), names = NA, las = las, lwd = 1, cex.lab=1.2, at=x, xaxt='n')
  for (i in 2:n) {
      boxplot(stats[[i]], outline=ifelse(plotStrip, FALSE, TRUE), col=cols[i],
              shrink=0.8, lwd=1, cex.lab=1.2, at=(x+i-1), add=TRUE, axes=FALSE)
  }
  
  if (plotGrid) {
    #grid(col=dark$gray)
    abline(h=axTicks(2), lty=3, col=mPal(oxCols$black))
  }
  
  if (plotStrip) {
      for (i in 1:n) {
          stripchart(stats[[i]], vertical = TRUE, method = "jitter", jitter = 0.2, add = TRUE, pch = 16, col = mPal(oxCols$red2,0.2), cex=0.8, at=(x+i-1))
      }
  }
  
  text(x=x+0.5, y=ylim[1]-0.25*(ylim[2]-ylim[1]), srt=45, labels=names, xpd=TRUE, cex=1.2)
}




plotPointProcess <- function(points, y1, y2, col=mPal(oxCols$blue1), ...) {
  
  n <- length(points)
  for (j in 1:n) {
    lines( rep(points[j], 2), c(y1, y2), col=col, ...)
  }
  points(points, rep(y2, n), pch=16, col=col, ...)
  
}



#' Plot a summary of a single simulation replicate's results
#' 
#' @param phylo
plotSimResult <- function(phylo, 
                          popSizeHPD, popSizeTimes, popSizeTraj, popSizeTruth,
                          sampIntensityHPD=NULL, sampIntensityTimes=NULL, sampIntensityTruth=NULL,
                          logPopSize=TRUE, short=TRUE, outputfile="") {
  
    if (outputfile != "") {
        pdf(outputfile, width=7.5, height=10)
    }
    layout(matrix(c(1,2,3), ncol=1), heights = c(2,1,1))

    k <- phylo$Nnode
    coalTimes <- beastio::getBranchingTimes(phylo)
    sampTimes <- beastio::getSamplingTimes(phylo)
    xlims <- c(0, max(coalTimes)) 
    
    
    
    # Plot tree
    par(mar=c(2,4,1,1))
    plot(ladderize(phylo, right=TRUE),  show.tip.label=FALSE, direction="leftwards", edge.width=0.5, edge.color=mPal(oxCols$gray6),
         xaxs='i', yaxs='i', xpd=TRUE, y.lim = c(-0.25*k, 1.04*k), x.lim=xlims)
    #nodelabels(pch=16, col=mPal(oxCols$orange1,0.5), cex=0.8, xpd=TRUE)
    #tiplabels(pch=16,  col=mPal(oxCols$blue1,0.5),   cex=0.8, xpd=TRUE)
    #rect(0, -0.25*k, max(coalTimes), 1.04*k, xpd=TRUE)
    
    
    
    # Plot sampling and coalescent events
    plotPointProcess(sampTimes, -0.1*k, -0.05*k, col=mPal(oxCols$blue1, 0.25), xpd=TRUE)
    axis(1, pos = -0.1*k, lwd = 0, lwd.ticks = 1, at=xlims, labels=NA)
    axis(1, pos = -0.1*k, at=c(axTicks(1), xlims[2]), labels=c(axTicks(1),NA))
    #axis(2, at = c(-0.1, -0.05)*k, labels=c(0,1), las=1)
    
    plotPointProcess(coalTimes, -0.25*k, -0.2*k, col=mPal(oxCols$gray3, 0.25), xpd=TRUE)
    axis(1, pos = -0.25*k, lwd = 0, lwd.ticks = 1, at=xlims, labels=NA)
    axis(1, pos = -0.25*k, at=c(axTicks(1), xlims[2]), labels=c(axTicks(1),NA))
    #axis(2, at = c(-0.25, -0.2)*k, labels=c(0,1), las=1)
    
    nrX <- xlims[1] - 0.06*diff(xlims)
    mtext("A", side=3, line=-1, at=nrX, cex=1.2)
  
    legend("bottomright", horiz=TRUE, inset=c(0, 0.15), legend=c("Sampling events", "Coalescent events"), col=c(mPal(oxCols$blue1), mPal(oxCols$gray3)), pch=16, bty='n')
    
    
    
    if (short) {
        #axis(1, pos = -0.25*k, lwd = 0, lwd.ticks = 1)
        #lines(xlims, rep(-0.25*k,2), xpd=TRUE)
      
        xlims <- c(0,max(sampTimes))
        abline(v=xlims[2], lty=2, col=mPal(oxCols$red2)) 
    }
    
    
    # Plot population size 
    if (logPopSize) {
        popSizeHPD   <- log(popSizeHPD)
        popSizeLabel <- "log N"
        popSizeTruth <- log(popSizeTruth)
        traj         <- function(x) log(popSizeTraj(x))
    } else {
        popSizeLabel <- "N"
        traj         <- popSizeTraj
    }  
    
    par(mar=c(2,4,2,1))
    plotSkyline(popSizeTimes, t(popSizeHPD), type='step', xlim=xlims, xaxs='i', yaxs='i', las=1, 
                xlab="", ylab=popSizeLabel, fill=mPal(oxCols$gray3, 0.5))
    x <- seq(xlims[1], xlims[2], by=0.01)
    
    # Plot true trajectory
    lines(x, traj(x), type="l", col=mPal(oxCols$green1), xaxs='i', lty=2, lwd=1)
    
    # Plot interval harmonic mean
    for (j in 2:length(popSizeTimes)) {
        lines(popSizeTimes[c(j-1, j, j)], popSizeTruth[c(j-1, j-1, j)], col=mPal(oxCols$red2), lty=1)
    }
    
    legend("topright", horiz=TRUE, inset=c(0,-0.12), legend=c("True N", "Segment harmonic mean"), col=c(mPal(oxCols$green1), mPal(oxCols$red2)), 
           lty=c(2,1), bty='n', xpd=TRUE)
    
    nrX <- xlims[1] - 0.06*diff(xlims)
    mtext("B", side=3, line=0, at=nrX, cex=1.2)
    
    
    
    # Plot sampling intensity
    if (!is.null(sampIntensityHPD)) {
        
        if (!is.matrix(sampIntensityHPD)) {
            sampIntensityHPD <- as.matrix(t(sampIntensityHPD))
        }
      
        if (any(is.na(sampIntensityTruth))) {
            ylims <- c(0, max(sampIntensityHPD))*1.05
        } else {
            ylims <- c(0, max(sampIntensityTruth, sampIntensityHPD))*1.05  
        }
        
        plotSkyline(sampIntensityTimes, t(sampIntensityHPD), type='step', xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', las=1, 
                    xlab="Time (from present)", ylab=expression(beta), fill=mPal(oxCols$gray3, 0.5))
        
        if (all(!is.na(sampIntensityTruth))) {
            for (j in 2:length(sampIntensityTimes)) {
                lines(sampIntensityTimes[c(j-1, j, j)], sampIntensityTruth[c(j-1, j-1, j)], col=mPal(oxCols$red2))
            }
            legend("topright", inset=c(0,-0.12), legend=expression("True"~beta), col=mPal(oxCols$red2), 
                   lty=1, bty='n', xpd=TRUE)
        }
        mtext("C", side=3, line=0, at=nrX, cex=1.2)
    }
    
    
    if (outputfile != "") {
        dev.off()
    }
  
}



