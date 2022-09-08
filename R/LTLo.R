#'
#'The Logarithmic Transformed Lomax
#'
#'Density, survival function, hazard rate function, distribution function, quantile
#'function and random generation for the Logarithmic Transformed Lomax distribution.
#'This distribution is an extension of the traditional Lomax distribution with one
#'scale parameter and two shape parameters.
#'
#'
#'@param \code{x,q} vector of quantiles.
#'@param \code{p} vector of probabilities.
#'@param \code{n} number of observations. If length(n) > 1, the length is taken
#'to be the number required.
#'@param \code{phi,sigma} vectors of shape parameters.
#'@param \code{tau} vector of scale parameters.
#'@param \code{log,log.p} logical; if TRUE, probabilities p are given as log(p).
#'@param \code{lower.tail} logical; if TRUE (default), probabilities are P(X≤x), otherwise, P(X>x).
#'
#'@details If tau is not specified, it assumes the default value of 1.
#'@details The Logarithmic Transformed Lomax distribution has density
#'@details \deqn{f(x;\phi,\sigma,\tau) = \frac{(\phi-1)\sigma\tau^\sigma(\tau+x)^
#'{-(\sigma+1)}}{\log(\phi)[1-(1-\phi)\tau^\sigma(\tau^\sigma(\tau+x)^{-\sigma})]}}
#'@details for \eqn{x,\phi,\sigma,\tau \geq 0  e \phi \neq 1.}
#'
#'@return dLTLo gives the density, sLTLo gives the survival function, hLTLo gives
#'the hazard rates function, pLTLo gives the distribution function, qLTLo gives the
#'quantile function, and rLTLo generates random deviates.
#'@return The length of the result is determined by n for rLTo, and is the maximum of
#'the lengths of the numerical arguments for the other functions.
#'@return The numerical arguments other than n are recycled to the length of the result.
#'Only the first elements of the logical arguments are used.
#'
#'@references Alotaibi, R., Okasha, H., Rezk, H., Almarashi, A. M., & Nassar, M. (2021).
#' On a new flexible Lomax distribution: statistical properties and estimation procedures
#' with applications to engineering and medical data. AIMS Mathematics, 6(12), 13976-13999.
#'
#'Lomax, K. S. (1954). Business failures: Another example of the analysis of failure data.
#'Journal of the American Statistical Association, 49(268), 847-852.
#'
#'@import stats
#'
#'@examples
#'
#'
# Densidade - Gráfico
#'x<-seq(0,2,length=100)
#'fx1 <- dLTLo(x=x,phi=20,sigma=20,tau=5)
#'fx2 <- dLTLo(x=x,phi=50,sigma=20,tau=5)
#'fx3 <- dLTLo(x=x,phi=2,sigma=5,tau=3)
#'fx4 <- dLTLo(x=x,phi=0.05,sigma=0.9,tau=2)
#'
#'plot(x,fx1,type="l",col=1,yaxp=c(0,2,10),ylab=c("f(x)"))
#'points(x,fx2,type="l",col=2)
#'points(x,fx3,type="l",col=3)
#'points(x,fx4,type="l",col=4)
#'legend(1.35,1.2,lty=1, col=c(1,2,3,4),lwd=1,c("φ=20     σ=20   τ=5","φ=50     σ=20   τ=5","φ=2       σ=5      τ=3","φ=0.05  σ=0.9  τ=2"))
#'
# taxa de falha - gráfico
#'x <- seq(0,3,length=100)
#'hx1 <- hLTLo(x=x, phi=50,  sigma=20)
#'hx2 <- hLTLo(x=x, phi=100, sigma=3.8)
#'hx3 <- hLTLo(x=x, phi=0.02,sigma=2)
#'
#'par(mfrow=c(1,3))
#'plot(x,hx1,type="l",col=1,xlim=c(0,3))
#'plot(x,hx2,type="l",col=1,xlim=c(0,3))
#'plot(x,hx3,type="l",col=1,xlim=c(0,3))
#'
#'
#'
#'
#'
#'@export
dLTLo <- function(x, phi, sigma, tau = 1, log = FALSE){
  stopifnot(phi>0 && sigma>0 && tau>0 && phi!=1)
  fun <- ((phi-1)*sigma*(tau^sigma)*(tau+x)^(-sigma-1))/
    (log(phi)*(1-(1-phi)*(tau^sigma)*(tau+x)^(-sigma)))
  if (log == TRUE){
    fun <- log(fun)
  }
  return(fun)
}
#' @rdname dLTLo
sLTLo <- function(x, phi, sigma, tau = 1, log = FALSE){
  stopifnot(phi>0 && sigma>0 && tau>0 && phi!=1)
  fun <- (log(phi)*(1-(1-phi)*(tau^sigma)*(tau+x)^(-sigma)))/log(phi)
  if (log == TRUE){
    fun <- log(fun)
  }
  return(fun)
}
#' @rdname dLTLo
hLTLo <-function(x, phi, sigma, tau = 1, log = FALSE){
  stopifnot(phi>0 && sigma>0 && tau>0 && phi!=1)
  fun <- ((phi-1)*sigma*(tau^sigma)*(tau+x)^(-sigma-1))/
    ((1-(1-phi)*(tau^sigma)*(tau+x)^(-sigma))*log(1-(1-phi)*(tau^sigma)*(tau+x)^(-sigma)))
  if (log == TRUE){
    fun <- log(fun)
  }
  return(fun)
}
#' @rdname dLTLo
pLTLo <-function(q, phi, sigma,tau, lower.tail = TRUE, log.p = FALSE){
  stopifnot(phi>0 && sigma>0 && tau>0 && phi!=1)
  fun <- 1-(log(1-(1-phi)*(tau^sigma)*(tau+q)^(-sigma))/log(phi))
  if (lower.tail == FALSE){
    fun <- 1-fun
  }
  if (log.p == TRUE){
    fun <- log(fun)
  }
  return(fun)
}
#' @rdname dLTLo
qLTLo <-function(q, phi, sigma, tau, lower.tail = TRUE, log.p = FALSE){
  stopifnot(phi>0 && sigma>0 && tau>0 && phi!=1)
  fun <- ((1-exp(log(phi)*(1-q)))/((1-phi)*tau^sigma))^(-1/sigma)-tau
  if (lower.tail == FALSE){
    fun <- 1-fun
  }
  if (log.p == TRUE){
    fun <- log(fun)
  }
  return(fun)
}
#' @rdname dLTLo
rLTLo <- function(n, phi, sigma, tau){
  stopifnot(phi>0 && sigma>0 && tau>0 && phi!=1)
  u   <- runif(n)
  fun <- ((1-exp(log(phi)*(1-u)))/((1-phi)*tau^sigma))^(-1/sigma)-tau
  return(fun)
}
