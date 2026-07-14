# KEEP

#' @title Adaptive Gaussian-Hermite variance estimation
#'
#' @description
#' Computes standard deviation estimates using adaptive Gaussian-Hermite
#' quadrature for a nonlinear mixed-effects model.
#'
#' @param RespLog A list containing symbolic log-likelihood components.
#' @param long.data A longitudinal dataset containing model observations.
#' @param idVar Character string specifying the subject identifier variable.
#' @param fixedest0 Estimated fixed-effect parameters.
#' @param dispest0 Estimated dispersion parameters.
#' @param invSIGMA0 Inverse of the random-effects covariance matrix.
#' @param Bi Estimated subject-specific random effects.
#' @param B Matrix of random-effect estimates.
#' @param Jfixed Indices or definitions of fixed-effect parameters.
#' @param Jraneff Indices or definitions of random-effect parameters.
#' @param ghsize Number of Gaussian-Hermite quadrature points. Default is \code{4}.
#' @param Silent Logical indicating whether errors are suppressed. Default is
#' \code{TRUE}.
#' @param epsilon Numerical tolerance for quadrature calculations. Default is
#' \code{10^{-6}}.
#' @param parallel Logical indicating whether parallel computation is used.
#' Default is \code{FALSE}.
#'
#' @return Standard deviation estimates obtained from adaptive Gaussian-Hermite
#' quadrature.
#'
#' @details
#' Generates subject-specific Gaussian-Hermite quadrature samples using estimated
#' random-effects covariance matrices and computes variance estimates through
#' numerical integration.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   # Example:
#'   # sd <- get_sd_aGH(...)
#' }
#' }
#'
#' @rdname get_sd_aGH
#' @export
get_sd_aGH <- function(RespLog, long.data, idVar, 
                       fixedest0, dispest0, invSIGMA0,Bi, B,
                       Jfixed, Jraneff,  
                       ghsize=4, Silent=T, epsilon=10^{-6}, 
                       parallel=F){
  q <- ncol(Bi)
  group <- long.data[ , idVar]  
  uniqueID <- unique(group)   
  
  GHzsamp0 = mgauss.hermite(n=ghsize, mu=rep(0,q), sigma=NULL)
  
  idSIGMA = get_idSIGMA_aGH(RespLog, long.data, idVar, uniqueID,
                            fixedest0, dispest0, invSIGMA0,Bi, B,
                            Jfixed, Jraneff) 
  # generate GH samples by subject
  GHsample0 <- as.list(rep(NA,n))
  for(i in 1:n){
    GHsample0[[i]] = mgauss.hermite(n=ghsize, mu=as.numeric(Bi[i,]), sigma=idSIGMA[[i]])
  }
  
  GHsd2 = try(calculate_aGH(RespLog, long.data, idVar, uniqueID,
                            fixedest0, dispest0, invSIGMA0,
                            GHzsamp0,GHsample0,
                            Jfixed, Jraneff,
                            ghsize, epsilon, parallel),  silent = Silent)
  return(GHsd2)
}
