#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param RespLog PARAM_DESCRIPTION
#' @param long.data PARAM_DESCRIPTION
#' @param idVar PARAM_DESCRIPTION
#' @param fixedest0 PARAM_DESCRIPTION
#' @param dispest0 PARAM_DESCRIPTION
#' @param invSIGMA0 PARAM_DESCRIPTION
#' @param Bi PARAM_DESCRIPTION
#' @param B PARAM_DESCRIPTION
#' @param Jfixed PARAM_DESCRIPTION
#' @param Jraneff PARAM_DESCRIPTION
#' @param ghsize PARAM_DESCRIPTION, Default: 4
#' @param Silent PARAM_DESCRIPTION, Default: T
#' @param epsilon PARAM_DESCRIPTION, Default: 10^{
#'    -6
#'}
#' @param parallel PARAM_DESCRIPTION, Default: F
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
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
