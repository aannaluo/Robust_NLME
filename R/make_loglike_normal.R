#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param par PARAM_DESCRIPTION
#' @param mean PARAM_DESCRIPTION
#' @param sd PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname make_loglike_normal
#' @export 
make_loglike_normal <- function(par, mean, sd){
  len <- length(par)
  
  temp_loglike <- as.list(rep(NA,len))
  
  for(i in 1:len){
    temp_par <- par[i]
    temp_mean <- mean[i]
    temp_sd <- sd[i]
    temp_loglike[[i]] <- paste0("-0.5*(", temp_par, "-", temp_mean, 
                                ")^2/",temp_sd,"^2", "-log(",temp_sd,")-0.5*log(2*pi)")
  }
  
  loglike <- paste0("(", unlist(temp_loglike), ")", collapse="+") # assume independent 
  return(loglike)
}
