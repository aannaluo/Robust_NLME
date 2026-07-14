#' @export
summary.Rnlme <- function(object, ...) {
  
  cat("Joint nonlinear mixed-effects model\n")
  cat("----------------------------------\n")
  
  cat("Convergence:",
      ifelse(object$convergence, "Successful", "Failed"),
      "\n\n")
  
  cat("Fixed effects:\n")
  print(object$fixedest)
  
  cat("\nDispersion parameters:\n")
  print(object$dispersion)
  
  invisible(object)
}

#' @export
coef.Rnlme <- function(object, ...) {
  object$fixedest
}

#' @export
AIC.Rnlme <- function(object, ...) {
  object$AIC
}

#' @export
BIC.Rnlme <- function(object, ...) {
  object$BIC
}

#' @export
var.Rnlme <- function(x, ...) {
  
  object <- x
  
  object$SIGMA
}

#' @export
ranef.Rnlme <- function(object, ...) {
  
  object$Bi
}