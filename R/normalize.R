#' normalize a series using conditional moments
#'
#' This function produces  a normalized series using conditional moments.
#' @param y The variable name
#' @param fit_mean 	Model object of class "conditional_moment" returned from
#'  \code{\link[conduits]{conditional_mean}} with information to append to observations.
#' @param fit_var 	Model object of class "conditional_moment" returned from
#'  \code{\link[conduits]{conditional_var}} with information to append to observations.
#' @param data a tibble containing all the time series
#' which are uniquely identified by the corresponding
#' Timestamp.
#' @return A \code{\link[tibble]{tibble}} with the conditional normliased series
#' @importFrom dplyr mutate bind_cols
#' @importFrom stats predict qnorm
#' @importFrom broom augment
#' @examples
#' data <- NEON_PRIN_5min_cleaned %>%
#' dplyr::filter(site == "upstream") %>%
#' dplyr::select(Timestamp, turbidity, level, conductance, temperature)
#'
#' fit_mean <- data %>%
#' conditional_mean(turbidity ~ s(level, k = 8 ) +
#' s(conductance, k = 8) + s(temperature, k = 8))
#'
#' fit_var <- data %>%
#' conditional_var(turbidity ~ level + conductance + temperature, fit_mean,
#' knots = c(7,7,7), family = "Gamma")
#'
#' #y_norm <-  data %>% normalize( turbidity, fit_mean, fit_var)
#'
#' @export
#'
normalize <- function(data, y, fit_mean, fit_var){

  y <- dplyr::enquo(y)
  cond_EY <- as.numeric(mgcv::predict.gam(fit_mean, newdata = data))
  cond_VY <- as.numeric(mgcv::predict.gam(fit_var, newdata = data,
                                           type = "response"))
  y_norm <- (data[y] - cond_EY)/ sqrt(cond_VY)
  return(y_norm)
}

