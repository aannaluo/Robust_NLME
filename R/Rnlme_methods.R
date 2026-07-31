#' @export
summary.Rnlme <- function(object, ...) {
  
  estimates <- object$fixedest
  se        <- object$fixedSD
  
  if (is.null(names(se)) && length(se) == length(estimates)) {
    names(se) <- names(estimates)
  }
  
  fixed_role <- sub("[0-9]+$", "", names(estimates))
  
  z <- as.numeric(estimates) / as.numeric(se)
  p <- 2 * pnorm(-abs(z))
  
  fixed_df <- data.frame(
    label = names(estimates),
    role  = fixed_role,
    est   = as.numeric(estimates),
    se    = as.numeric(se),
    `z value`  = z,
    `p-value` = p,
    check.names = FALSE,
    row.names = names(estimates)
  )
  
  # ---- dispersion parameters (d, siga, sigb, xi, ...) ----
  disp_est <- unlist(object$dispersion)
  disp_se  <- object$dispSD
  if (is.null(names(disp_se)) && length(disp_se) == length(disp_est)) {
    names(disp_se) <- names(disp_est)
  }
  disp_role <- sub("[0-9]+$", "", names(disp_est))
  
  disp_df <- NULL
  if (length(disp_se) == length(disp_est)) {
    dz <- as.numeric(disp_est) / as.numeric(disp_se)
    dp <- 2 * pnorm(-abs(dz))
    disp_df <- data.frame(
      label = names(disp_est),
      role  = disp_role,
      est   = as.numeric(disp_est),
      se    = as.numeric(disp_se),
      `z value`  = dz,
      `p-value` = dp,
      check.names = FALSE,
      row.names = names(disp_est)
    )
  } else {
    disp_df <- data.frame(
      label = names(disp_est),
      role  = disp_role,
      est   = as.numeric(disp_est),
      check.names = FALSE,
      row.names = names(disp_est)
    )
  }
  
  # ---- formulas per model component (unchanged) ----
  formulas <- lapply(object$nlmeObjects, function(obj) {
    list(
      mean_model = obj$model,
      fixed      = obj$fixed,
      random     = obj$random,
      disp_model = if (!is.null(obj$sigma)) obj$sigma$model else NULL
    )
  })
  
  out <- list(
    convergence = object$convergence,
    formulas    = formulas,
    fixed_df    = fixed_df,
    disp_df     = disp_df,
    AIC         = object$AIC,
    BIC         = object$BIC,
    logLik      = object$loglike_value,
    corr        = object$SIGMA
    )
  
  class(out) <- "summary.Rnlme"
  out
}
  
  # p-value
  # use a normal distribution
  # just for the fixed effects --> 
  # 
  # standard error
  
  # estimate <- 0.75
  # se <- 0.20
  
  # z <- estimate / se
  # p <- 2 * pnorm(-abs(z))
  
  # z
  # p
  
  # separate using dashed lines, mean model (beta), variance model (alpha), 
  # measurement error model (gamma) in t-table
  
  # simple formula for each model
  
  

#' @export
print.summary.Rnlme <- function(x, digits = 4, ...) {
  
  cat("Joint Nonlinear Mixed-Effects Model\n")
  cat(sprintf("logLik: %.2f   AIC: %.2f   BIC: %.2f\n",
              x$logLik, x$AIC, x$BIC))
  cat("-----------------------------------\n\n")
  
  print_block <- function(header, formula_lines, df, role) {
    cat(header, "\n")
    for (nm in names(formula_lines)) {
      if (!is.null(formula_lines[[nm]])) {
        cat("  ", nm, ": ", deparse(formula_lines[[nm]]), "\n", sep = "")
      }
    }
    cat("\n")
    sub_df <- df[df$role == role, , drop = FALSE]
    if (nrow(sub_df) > 0) {
      mat <- as.matrix(sub_df[, c("est", "se", "p-value")])
      rownames(mat) <- sub_df$label
      print(mat, digits = digits)
    }
    cat("\n-----------------------------------\n\n")
  }
  
  f1 <- x$formulas[[1]]
  print_block("Mean model (beta):",
              list(model = f1$mean_model, fixed = f1$fixed, random = f1$random),
              x$fixed_df, "beta")
  
  print_block("Variance model (alpha):",
              list(model = f1$disp_model),
              x$fixed_df, "alpha")
  
  if (length(x$formulas) > 1) {
    f2 <- x$formulas[[2]]
    print_block("Measurement error model (gamma):",
                list(model = f2$mean_model, fixed = f2$fixed, random = f2$random),
                x$fixed_df, "gamma")
  }
  
  cat("Dispersion parameters:\n\n")
  disp_mat <- as.matrix(x$disp_df[, c("est", "se", "p-value")])
  rownames(disp_mat) <- x$disp_df$label
  print(disp_mat, digits = digits)
  cat("\n")
  
  
  invisible(x)
}

#' @export
fixef.Rnlme <- function(object, ...) {
  object$fixedest
  
  #output t-table for coefficients from summary function.
  # subset of the summary function
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
VarCorr.Rnlme <- function(object, ...) {
  
  SIGMA <- object$SIGMA
  
  ran_names <- unlist(lapply(object$nlmeObjects, function(obj) {
    rnames <- all.vars(obj$random[[2]])
    paste0(obj$ranName, ".", rnames)
  }))
  
  if (length(ran_names) == nrow(SIGMA)) {
    rownames(SIGMA) <- ran_names
    colnames(SIGMA) <- ran_names
  }
  
  SIGMA
  
  # add labels.
}

#' @export
ranef.Rnlme <- function(object, ...) {
  
  Bi <- object$Bi
  
  out <- data.frame(
    id = object$uniqueID,
    Bi,
    check.names = FALSE
  )
  
  rownames(out) <- NULL
  out
  
  # add id column
}

#' @export
loglik.Rnlme <- function(object, ...) {
  object$loglike_value
}

#' @export 
fitted.Rnlme <- function(object, ...) {
  # estimates of fixed effects and random effects
  # calculate values for patient 
  
  
}
