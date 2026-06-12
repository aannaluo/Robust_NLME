# evaluate a matrix 

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mat PARAM_DESCRIPTION
#' @param dim PARAM_DESCRIPTION
#' @param data PARAM_DESCRIPTION, Default: NULL
#' @param par.val PARAM_DESCRIPTION, Default: NULL
#' @param raneff.val PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname get_Hvalue
#' @export 
get_Hvalue <- function(mat, dim, data=NULL, par.val=NULL, raneff.val=NULL){
  D <- matrix(0, dim, dim)
  for(i in 1:dim){
    for(j in 1:dim){
      kk <- (i-1)*dim+j
      Di <- with(data, with(par.val, with(raneff.val,  
                                          eval(parse(text=mat[[kk]])))))
      
      if(length(Di)==1 & !is.null(raneff.val)) Di <- Di*nrow(raneff.val)
      D[i,j] <- sum(Di, na.rm=TRUE)
    }
  }
  
  return(D)
}


