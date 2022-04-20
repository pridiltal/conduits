globalVariables(c(".fitted", ".se.fit", ".cond_m"))

#' Augment data with information from an conditional mean fit
#'
#' This function produces partial residuals for each predictor,
#' and the estimated conditional means, standard error and confidence limits.
#'
#' @param x 	Model object of class "conditional_moment" returned from
#'  \code{\link[conduits]{conditional_mean}} or  \code{\link[conduits]{conditional_var}}
#'  with information to append to observations.
#' @param level 	Confidence level. Default is set to 0.95.
#' @param ...	 Addition arguments to augment method.
#' @return A \code{\link[tibble]{tibble}} with information
#'  about data points.
#'
#' @seealso \code{\link[mgcv]{gam}}
#' @importFrom dplyr mutate bind_cols rename
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
#' data_inf <- fit_mean %>% augment()
#'
#' @export augment.conditional_moment
#' @export
#'
augment.conditional_moment <- function(x, level = 0.95, ...){

  # getting each component of the linear predictor from the fitted model
  fv <- stats::predict(x, type = "terms")

  aug <- broom:::augment.gam(x) # check before cran

  .partial.res = fv + aug$.resid
  colnames(.partial.res) <- paste0('.presid_', colnames(.partial.res))

  zq <- abs(stats::qnorm((1-level)/2))
  data <- aug %>%
    dplyr::mutate(
      .cond_m = .fitted,
      .LI = as.numeric(.fitted - zq*.se.fit),
      .UI = as.numeric(.fitted + zq*.se.fit),
    ) %>%
    dplyr::bind_cols(.partial.res)


  #if(inherits(x, "conditional_mean"))
  if(x$type == "conditional_mean")
    data <- data %>% dplyr::rename(.cond_EX = .cond_m)
  #if(inherits(x, "conditional_var"))
  if(x$type == "conditional_var")
    data <- data %>% dplyr::rename(.cond_VAR = .cond_m)

  return(data)
}


#' @export
broom::augment
