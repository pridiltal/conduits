#' Estimating conditional variance of a time series
#'
#' This function estimates the variance of a time
#' series conditional on a set of  other times series
#' via additive models.
#'
#' @param data A tibble containing all the time series
#' which are uniquely identified by the corresponding
#' Timestamp.
#' @param formula An object of class "formula": a symbolic description of the model to be fitted.
#' The details of model specification are given under ‘Details’.
#' @param fit_mean A GAM object return from  \code{\link[conduits]{conditional_mean}}
#' @param family the family to be used in conditional variance model. Currently
#' this can take either "Gamma" or "lognormal".
#' @param knots a vector specifying the dimension of the basis in the smooth term fitting for
#' each predictor in the GAM for conditional variance of $x$. Each component of the vector should
#' corresponds to each numeric predictor specified in the formula
#' @return The function returns an object of class
#' "gam" as described in \code{\link[mgcv]{gamObject}}.
#' @details{ Suppose $x_t$ is a time series where its
#' variance is a function of $z_t$. i.e. $Var(x_t|z_t) = v_x(z_t)$.
#' Then $v_x(z_t)$can be estimated via generalised
#' additive models (GAM). This function uses GAMs implemented
#' in \code{mgcv} package to estimate the conditional variance
#' of a time series given a set of time series predictors.}
#'
#' @seealso \code{\link[mgcv]{gam}}
#' @importFrom mgcv predict.gam
#' @importFrom mgcv gam
#' @importFrom stats as.formula
#' @importFrom stats Gamma
#' @export
#' @examples
#'
#' data <- NEON_PRIN_5min_cleaned %>%
#' dplyr::filter(site == "upstream") %>%
#' dplyr::select(Timestamp, turbidity, level, conductance, temperature)
#'
#' fit_mean <- data %>%
#' conditional_mean(turbidity ~ s(level, k = 8 ) +
#' s(conductance, k = 8) + s(temperature, k = 8))
#'
#' \dontrun{
#' fit_var <- data %>%
#' conditional_var(turbidity ~ level + conductance + temperature, fit_mean,
#' knots = c(7,7,7), family = "Gamma")
#' }
conditional_var <- function(data, formula, fit_mean,
                            family, knots = NULL)
{
  vars <- all.vars(formula)
  z_fac <- data %>% Filter(f = is.factor) %>% names
  y_name <-  vars[1]
  y <- data[[y_name]]
  z_num <- vars[!(vars %in% c(y_name, z_fac))]

  # if variance knots are null replace with the default in s()
  if(is.null(knots)){
    knots_variance <- rep(-1, length(z_num))
  }

  if(rlang::is_empty(z_fac)){

    if(family == "Gamma"){
      formula_var <- paste("Y_Ey2 ~", paste("s(", z_num, ", k=", knots,  ")",
                                              sep = "", collapse = " + "), sep = " ")
    }
    if(family == "lognormal"){
      formula_var <- paste("log(Y_Ey2) ~", paste("s(", z_num, ", k=", knots,  ")",
                                                   sep = "", collapse = " + "), sep = " ")
    }

  } else{
    if(family == "Gamma"){
      formula_var <- paste("Y_Ey2 ~", paste("s(", z_num, ", k=", knots,  ")",
                                              sep = "", collapse = " + "),
                             "+", paste(z_fac, collapse = " + "),
                             sep = " ")
    }

    if(family == "lognormal"){
      formula_var <- paste("log(Y_Ey2) ~", paste("s(", z_num, ", k=", knots,  ")",
                                                   sep = "", collapse = " + "),
                             "+", paste(z_fac, collapse = " + "),
                             sep = " ")
    }

  }

  # Compute conditional means, squared errors from the x_mean_gam
  data <- data %>%
    dplyr::mutate(E_Y = as.numeric(mgcv::predict.gam(fit_mean,
                                                     newdata = data)),
                  Y_Ey2 = (y - .data$E_Y)^2)

  if(family == "Gamma"){
    var_gam_fit <- mgcv::gam(formula = stats::as.formula(formula_var),
                           data = data,
                           family = stats::Gamma(link = "log"))
  }

  if(family == "lognormal"){
    var_gam_fit <- mgcv::gam(formula = stats::as.formula(formula_var),
                           data = data,
                           family = stats::gaussian())
  }

  class(var_gam_fit) <- c("conditional_moment", "conditional_var", "gam", "glm", "lm" )
  return(var_gam_fit)
}
