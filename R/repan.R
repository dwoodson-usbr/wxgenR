#' Random variates from the Epanechnikov kernel
#'
#' Simulate outside the historical envelope
#' using randomly generated values from the Epanechnikov kernel
#' (via acceptance-rejection sampling). \cr
#' \cr
#' For more details on the Epanechnikov kernel and its use
#' in a weather generator, see Rajagopalan et al. (1996).
#'
#' @param nsim Number of simulations.
#'
#' @return Returns a vector of random variates sampled from the Epanechnikov kernel. `nsim` number of samples are returned.
#'
#'
#' @examples
#'  repan(nsim = 10)
#'
#'  #simulate and plot density and distribution function
#'  oldpar = par(mfrow=c(1,3), mar=c(2,2.5,2,1),
#'               oma=c(2,2,0,0), mgp=c(2,1,0), cex.axis=0.8)
#'
#'  par(mfrow=c(1,2))
#'  nsim=1e5
#'  x <- sort(repan(nsim));y=0.75*(1-x^2)
#'  plot(x,y,xlab="x",ylab="f(x)",type="l",lwd=2)
#'  grid()
#'  title (main="Epanechnikov PDF",cex.main=0.8)
#'  F=rank(x)/(nsim+1)
#'  plot(x,F,ylab="F(x)",type="l",lwd=2)
#'  grid()
#'  title (main="Epanechnikov CDF",cex.main=0.8)
#'
#'  dev.off()
#'
#'  par(oldpar)
#'
#' @references {Rajagopalan, B., Lall, U., & Tarboton, D. G. (1996). Nonhomogeneous Markov Model for Daily Precipitation. Journal of Hydrologic Engineering, 1(1), 33–40. https://doi.org/10.1061/(ASCE)1084-0699(1996)1:1(33)}
#'
#' @export
#'
# @rawNamespace import(stats, except = filter)
#'
#'
#'

"repan" <- function(nsim){
  #simulations using Epanechnikov kernel
  #using acceptance-rejection sampling
  icount=0
  x <- rep(NA,nsim)
  while (icount < nsim){
    u1=runif(1,-1,1)
    u2=runif(1,0,1)
    if (((3*u1^2+4*u2)<=3)){
      icount=icount+1
      x[icount]=u1
    }
  } #end while
  return(x)
} #end function


