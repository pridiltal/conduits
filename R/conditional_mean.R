#' Estimating conditional mean of a time series
#'
#' This function estimates the means  of a time
#' series conditional on a set of  other times series
#' via additive models.
#'
#' @param data a tibble containing all the time series
#' which are uniquely identified by the corresponding
#' Timestamp.
#' @param formula an object of class "formula", a
#' symbolic description of the model to be fitted.
#' The details of model specification are given
#' under ‘Details’.
#' @param knots_mean  a vector specifying the dimension
#' of the basis in the smooth term fitting for each
#' predictor in the GAM for conditional mean of $x$.
#'  Each component of the vector should corresponds
#'  to each predictor specified in the formula.
#' @return The function returns an object of class
#' "gam" as described in \code{\link[mgcv]{gamObject}}.
#' @details{ Suppose $x_t$ is a time series where its
#' mean is a function of $z_t$. i.e. $E(x_t|z_t) = m_x(z_t)$.
#' Then $m_x(z_t)$ can be estimated via generalised
#' additive models (GAM). This function uses
#' GAMs implemented in \code{mgcv} package to estimate
#' the conditional means of a time series given a set of
#' time series predictors.}
#'
#' @seealso \code{\link[mgcv]{gam}}
#' @importFrom rlang is_empty
#' @importFrom stats as.formula
#' @export
#' @examples
#' mean_fit <- NEON_PRIN_5min_cleaned %>%
#' dplyr::filter(site == "upstream") %>%
#' dplyr::select(turbidity, level, conductance, temperature) %>%
#' conditional_mean(turbidity ~ level + conductance + temperature,
#'                  knots_mean = c(8, 8, 8))
#'
conditional_mean <- function(data, formula, knots_mean = NULL)
{
  vars <- all.vars(formula)
  z_fac <- data %>% Filter(f = is.factor) %>% names
  y <-  vars[1]
  z_num <- vars[!(vars %in% c(y, z_fac))]

  # if mean knots are null replace with the default in s()
  if(is.null(knots_mean)){
     knots_mean <- rep(-1, length(z_num))
  }

  if(rlang::is_empty(z_fac))
  {
    formula_new <- paste(y, "~", paste("s(", z_num, ", k=", knots_mean,  ")",
                                               sep="", collapse="+"),sep="")

  } else{
    formula_new <- paste(y, "~", paste("s(", z_num, ", k=", knots_mean, ")", sep="", collapse="+"),
                            "+", paste(z_fac, collapse = " + "),
                            sep = " ")
  }

  mean_gam_fit <- mgcv::gam(
    formula =  stats::as.formula(formula_new),
    data = data)

  return(mean_gam_fit)
}
