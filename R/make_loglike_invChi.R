# this function return the log density for inverse-Chi distribution 
# par: vector of variables, p>=1
# df
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param par PARAM_DESCRIPTION
#' @param df PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname make_loglike_invChi
#' @export 
make_loglike_invChi <- function(par, df){
  len <- length(par)
  
  temp_loglike <- as.list(rep(NA,len))
  
  for(i in 1:len){
    temp_par <- par[i]
    temp_loglike[[i]] <- paste0("log((0.5*", df, ")^(0.5*", df, ")/gamma(0.5*",df,"))-0.5*",df,"*", temp_par,"-0.5*", df, "/exp(", temp_par, ")")
  }
  
  loglike <- paste0("(", unlist(temp_loglike), ")", collapse="+")
  return(loglike)
}

