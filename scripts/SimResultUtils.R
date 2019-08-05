
approxHarmonicMean <- function(traj, x1, x2, steps=100) {
  x <- seq(x1, x2, length.out=steps)
  return ( 1/mean( 1/traj(x) ) )
}

approxArithmeticMean <- function(traj, x1, x2, steps=100) {
  x <- seq(x1, x2, length.out=steps)
  return ( 1/mean( 1/traj(x) ) )
}

approxGeomCombMean <- function(traj, x1, x2, a=1, b=1, steps=100) {
  x <- seq(x1, x2, length.out=steps)
  return(sqrt( (a/mean( 1/traj(x)) * b*mean(traj(x)))))
}

approxRange <- function(traj, x1, x2, steps=100) {
  x <- seq(x1, x2, length.out=steps)
  res <- range(traj(x))
  names(res) <- c("min", "max")
  return(res)
}

# Mean of parameter size in each interval
getMeanTraj <- function(t, traj, meanFn, ...) {
  
  parTruth <- c()
  parRange <- c()
  for (j in 2:length(t)) {
    parTruth <- c(parTruth, meanFn(traj, t[j-1], t[j], ...))
    parRange <- rbind(parRange, approxRange(traj, t[j-1], t[j]))
  }
  
  return(list(truth=parTruth, range=parRange))
}

piecewiseFn <- function(x, vals, breaks) {
  if (length(x) > 1) {
    return(sapply(x, piecewiseFn, simresult[[1]]$samp_intensity, simresult[[1]]$epoch_times))
  } else {
    return( vals[max(1, which.min(breaks < x) - 1)] )
  }
}

getPiecewiseFn <- function(vals, breaks) {
  return (function(x) piecewiseFn(x, vals, breaks))
}


#' Get summary statistics for a single run (simulated tree and inferred mcmc chain)
#' 
#' @param t The times of interval boundaries. Should start with 0 and be of dimension (n+1)
#' @param truthmean The mean (or harmonic mean) of the true trajectory for each interval. Should be of dimension n
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
  
  # Mean relative absolute deviation adn mean relative credible interval width
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
getStats <- function(simresult, popSizes, sampIntensities, popSizeTimes, sampIntensityTimes, cutPopSizes=TRUE, plotResults=TRUE, outputpath="") {
  
  if (outputpath != "" && plotResults == TRUE) {
    dir.create(paste0(outputpath,"figures/"), recursive=TRUE, showWarnings=FALSE)
    dir.create(paste0(outputpath,"figures-log/"), recursive=TRUE, showWarnings=FALSE)
    dir.create(paste0(outputpath,"figures-short/"), recursive=TRUE, showWarnings=FALSE)
    dir.create(paste0(outputpath,"figures-short-log/"), recursive=TRUE, showWarnings=FALSE)
    
    outnames <- sapply(strsplit(chanames(popSizes), "/"), function(x) x[length(x)])
    
    outnames <- sapply(strsplit(outnames, "\\."), function(x) paste(x[1:(length(x)-1)], collapse="."))
  }
  
  popSizeStats <- sampIntensityStats <- c()
  for (i in 1:length(popSizes)) {
    
    # Estimated HPDs
    popSizeHPD       <- getHPDMedian(popSizes[[i]])
    sampIntensityHPD <- getHPDMedian(sampIntensities[[i]])
    
    # Populations size group times    
    t <- c(0, popSizeTimes[[i]])
    s <- c(0, sampIntensityTimes[[i]])
    
    
    if (cutPopSizes == TRUE && s[length(s)] < t[length(t)]) {
      j <- which.min(t < s[length(s)])
      t <- t[1:j]
      t[j] <- s[length(s)]
      popSizeHPD <- popSizeHPD[1:(j-1),]
    }

    # True trajectory
    traj    <- function(x) simresult[[i]]$traj(x + simresult[[i]]$timeshift)
    logtraj <- function(x) log( simresult[[i]]$traj(x + simresult[[i]]$timeshift) )
    
    # Harmonic mean of population size in each interval
    popSizeTruth <- c()
    popSizeRange <- c()
    for (j in 2:length(t)) {
      popSizeTruth <- c(popSizeTruth, approxHarmonicMean(traj, t[j-1], t[j]))
      popSizeRange <- rbind(popSizeRange, approxRange(traj, t[j-1], t[j]))
    }
    
    # True samping intensity in each interval
    sampIntensityTruth <- simresult[[i]]$samp_intensity
    sampIntensityRange <- cbind(sampIntensityTruth, sampIntensityTruth)
    colnames(sampIntensityRange) <- c("min","max")

    if (outputpath != "" && plotResults == TRUE) {
        plotSimResult(simresult[[i]]$phylo, simresult[[i]]$coal_times, simresult[[i]]$samp_times, 
                      popSizeHPD, t, traj, popSizeTruth, 
                      sampIntensityHPD, s, sampIntensityTruth, 
                      logPopSize=FALSE, short=FALSE, outputfile=paste0(outputpath, "figures/", outnames[i],".pdf"))
        
        plotSimResult(simresult[[i]]$phylo, simresult[[i]]$coal_times, simresult[[i]]$samp_times, 
                      popSizeHPD, t, traj, popSizeTruth, 
                      sampIntensityHPD, s, sampIntensityTruth, 
                      logPopSize=FALSE, short=TRUE, outputfile=paste0(outputpath, "figures-short/", outnames[i],".pdf")) 
        
        plotSimResult(simresult[[i]]$phylo, simresult[[i]]$coal_times, simresult[[i]]$samp_times, 
                      popSizeHPD, t, traj, popSizeTruth, 
                      sampIntensityHPD, s, sampIntensityTruth, 
                      logPopSize=TRUE, short=FALSE, outputfile=paste0(outputpath, "figures-log/", outnames[i],".pdf"))
        
        plotSimResult(simresult[[i]]$phylo, simresult[[i]]$coal_times, simresult[[i]]$samp_times, 
                      popSizeHPD, t, traj, popSizeTruth, 
                      sampIntensityHPD, s, sampIntensityTruth, 
                      logPopSize=TRUE, short=TRUE, outputfile=paste0(outputpath, "figures-short-log/", outnames[i],".pdf")) 
    }
    
    # Calculate statistics  
    popSizeStats       <- rbind(popSizeStats,       getRunStats(t, popSizeTruth, popSizeRange, popSizeHPD))
    #sampIntensityStats <- rbind(sampIntensityStats, getRunStats(0:length(sampIntensityTruth), sampIntensityTruth, sampIntensityRange, sampIntensityHPD))
    sampIntensityStats <- rbind(sampIntensityStats, getRunStats(s, sampIntensityTruth, sampIntensityRange, sampIntensityHPD))
  }
  
  if (outputpath != "") {
      if (cutPopSizes == TRUE) {
          write.csv(popSizeStats,       paste0(outputpath, "popSizeStatsSamplingPeriod.csv"), row.names=FALSE)
      } else {
          write.csv(popSizeStats,       paste0(outputpath, "popSizeStatsFull.csv"), row.names=FALSE)
      }
      write.csv(sampIntensityStats, paste0(outputpath, "sampIntensityStats.csv"), row.names=FALSE)
  }
  
  return( list(popSizeStats = popSizeStats, 
               sampIntensityStats = sampIntensityStats))
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


plotPointProcess <- function(points, y1, y2, col=mPal(oxCols$blue1), ...) {
  
  n <- length(points)
  for (j in 1:n) {
    lines( rep(points[j], 2), c(y1, y2), col=col, ...)
  }
  points(points, rep(y2, n), pch=16, col=col, ...)
  
}



plotSimResult <- function(phylo, coalTimes, sampTimes, 
                          popSizeHPD, popSizeTimes, popSizeTraj, popSizeTruth,
                          sampIntensityHPD, sampIntensityTimes, sampIntensityTruth,
                          logPopSize=TRUE, short=TRUE, outputfile="") {
  
  if (outputfile != "") {
      pdf(outputfile, width=7.5, height=10)
  }
  xlims <- c(0, max(popSizeTimes)) 
  layout(matrix(c(1,2,3), ncol=1), heights = c(2,1,1))
  
  # Plot tree
  k <- phylo$Nnode
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
  
  
  #plotPointProcess(coalTimes, -0.15*k, -0.2*k, col=mPal(oxCols$gray3, 0.5), xpd=TRUE)
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
  
  # Interval harmonic mean
  for (j in 2:length(popSizeTimes)) {
      lines(popSizeTimes[c(j-1, j, j)], popSizeTruth[c(j-1, j-1, j)], col=mPal(oxCols$red2), lty=1)
  }
  
  # True trajectory
  lines(x, traj(x), type="l", col=mPal(oxCols$green1), xaxs='i', lty=2, lwd=1)
  
  legend("topright", horiz=TRUE, inset=c(0,-0.12), legend=c("True N", "Segment harmonic mean"), col=c(mPal(oxCols$green1), mPal(oxCols$red2)), 
         lty=c(2,1), bty='n', xpd=TRUE)
  
  
  nrX <- xlims[1] - 0.06*diff(xlims)
  mtext("B", side=3, line=0, at=nrX, cex=1.2)
  
  # Plot sampling intensity
  
  if (!is.matrix(sampIntensityHPD)) {
    sampIntensityHPD <- as.matrix(t(sampIntensityHPD))
  }
  
  plotSkyline(sampIntensityTimes, t(sampIntensityHPD), type='step', xlim=xlims, xaxs='i', yaxs='i', las=1, 
              xlab="Time (from present)", ylab=expression(beta), fill=mPal(oxCols$gray3, 0.5))
  for (j in 2:length(sampIntensityTimes)) {
      lines(sampIntensityTimes[c(j-1, j, j)], sampIntensityTruth[c(j-1, j-1, j)], col=mPal(oxCols$red2))
  }
  legend("topright", inset=c(0,-0.12), legend=expression("True"~beta), col=mPal(oxCols$red2), 
         lty=1, bty='n', xpd=TRUE)
  
  mtext("C", side=3, line=0, at=nrX, cex=1.2)
  
  
  if (outputfile != "") {
    dev.off()
  }
  
}