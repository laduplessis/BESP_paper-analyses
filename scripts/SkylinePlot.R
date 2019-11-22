###############################################################################
# Functions for plotting birth-death skylinesdev


#' Get range and pad by some percentage of the range on either side
#' Good for plot limits
paddedrange <- function(y, pad=0.1) {
  yrange <- c(min(y), max(y))
  return( c(yrange[1]-diff(yrange)*pad*0.5, yrange[2]+diff(yrange)*pad*0.5) )
}


#' Plot a Skyline.
#' 
#' Assume each column is a different time point (ncol(skylinemat) == length(times) must be fulfilled)
#' 
#' @param times The time points to draw the skyline at. length(times) should be equal to ncol(skyline_mat) + 1.
#'              (Because there are length(times) intervals and the last interval needs an end-time).
#'              If length(times) == ncol(skyline_mat) the last interval is arbitrarily given the same length as 
#'              the second last interval.
#' 
#' @param skyline_mat
#' 
#' @param type Type of skyline to plot.
#' 
#'             "smooth" : If skyline_mat contains the HPDs plot a smooth line and polygon (type='l')
#'             
#'             "step"   : If skyline_mat contains the HPDs plot a stepped line and polygon (type='S')
#'             
#'             "lines" : Plot every row in the matrix as a smooth line (type='l')
#'             
#'             "steplines" : Plot every row in the matrix as a stepped line (type='S')
#'             
#' @param traces Number of traces to draw if type="traces"
#'             
#' @param new Create a new set of axes (the skyline is fitted to the plotting devie)
#' 
#' @param add Add to the current plot (do not create a new plotting device)
#' 
#' @param ... Parameters passed to plotting function
#' 
#' @export
plotSkyline <- function(times, skyline_mat, type="smooth", traces=1000, col=mPal(oxCols$black), fill=mPal(oxCols$gray3, 0.25), lwd=1, lty=1,
                        new=TRUE, add=FALSE, xlims=NULL, ylims=NULL, ...) {
  
  ##################
  # Initialization #
  ##################
  
  intervals <- ncol(skyline_mat)
  skyline_mat <- cbind(skyline_mat, skyline_mat[,intervals])
  
  # Check dimensions
  if (length(times) == intervals)
      times <- c(times, times[intervals]+(times[intervals]-times[intervals-1]))
  else
  if (length(times) < intervals || length(times) > (intervals+1))
      stop("Dimension mismatch between times and skyline_mat!")
  
  # Do not open a new device (add to the current plot)
  if (add == TRUE) par(new=TRUE)
  
  # Create a new set of axes that fits this skyline
  if (new == TRUE) {
      # Always slightly pad the y-limits (even if provided) to plot values at the boundary properly
      if (is.null(ylims)) 
          ylims <- paddedrange(skyline_mat)
      else
          ylims <- paddedrange(ylims,0.01)
      if (is.null(xlims)) xlims <- range(times)
        
      # Create new plot
      plot(1, type='n', ylim=ylims, xlim=xlims, ...)
      
      # Get and set clipping area
      usr <- par("usr")
      clip(xlims[1], xlims[2], ylims[1], ylims[2])
  }
  
  #######################
  # Plot actual skyline #
  #######################
  
  # Plot traces
  if (type == "lines" || type == "steplines") {
      if (traces > 0 && traces < nrow(skyline_mat)) 
          ind <- sample(nrow(skyline_mat), size=traces, replace=FALSE)
      else 
          ind <- 1:nrow(skyline_mat)
    
      for (i in 1:length(ind)) {
        lines(times, skyline_mat[ind[i],], col=col, lwd=lwd, lty=lty, type=if(type=="lines") 'l' else 's')
      }   
  } else
  # Plot HPDs smooth
  if (type == "smooth" && nrow(skyline_mat) == 3) {
    
      polygon(c(times, rev(times)), c(skyline_mat[1,], rev(skyline_mat[3,])), col=fill, border=NA)
      lines(times, skyline_mat[2,], col=col, lwd=lwd, lty=lty)  
      
  } else
  # Plot HPDs stepped
  if (type == "step" && nrow(skyline_mat) == 3) {
    
      for (i in 2:(intervals+1)) {
        rect(times[i-1], skyline_mat[1,i-1], times[i], skyline_mat[3,i-1], col=fill, border=NA)
      }
      lines(times, skyline_mat[2,], col=col, lwd=lwd, lty=lty, type='s')
      
  } else
    stop("Invalid type parameter for Skyline plot!")
  
  ##################
  # Reset settings #
  ##################
  
  if (add == TRUE) par(new=FALSE)
  if (new == TRUE) do.call("clip", as.list(usr))
}




