#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param RespLog PARAM_DESCRIPTION
#' @param long.data PARAM_DESCRIPTION
#' @param fixedest PARAM_DESCRIPTION
#' @param dispest PARAM_DESCRIPTION
#' @param invSIGMA PARAM_DESCRIPTION
#' @param Bi PARAM_DESCRIPTION
#' @param B PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname get_loglike_value
#' @export 
get_loglike_value <- function(RespLog, long.data, fixedest, dispest, invSIGMA, Bi, B){
  
  Yrandisp <- !is.null(RespLog$randisp.loglike)
  Ysigma <- !is.null(RespLog$sigma.loglike)
  
  
  par.val <- as.list(c(fixedest, dispest))
  par.val$invSIGMA <- invSIGMA
  
  mu.val <- unlist(map(RespLog$mu.loglike, function(t){
    with(long.data, with(par.val, with(B,eval(parse(text=t)))))
  }))
  
  if(Ysigma){
    sigma.val <- with(par.val, with(Bi,eval(parse(text=RespLog$sigma.loglike))))
  } else sigma.val <- 0
  
  if(Yrandisp) {
    randisp.val <- with(par.val, with(Bi,eval(parse(text=RespLog$randisp.loglike))))
  } else {
    randisp.val <- 0
  }
  
  n <- nrow(Bi)
  ran.val <- vector("list", n)
  for(i in 1:n){
    ran.val[[i]] <-  with(par.val, with(Bi[i,], eval(parse(text=RespLog$ran.loglike))))
  }
  ran.val <- unlist(ran.val)
  
  return(sum(mu.val, na.rm=TRUE)+sum(sigma.val)+sum(randisp.val)+sum(ran.val))
}
