#' Augment data with information from an conditional mean fit
#'
#' This function produces partial residuals for each predictor,
#' and the estimated conditional means, standard error and confidence limits.
#'
#' @param x 	Model object returned from conditional_mean
#'  with information to append to observations.
#' @param ...	 Addition arguments to augment method.
#' @return A \code{\link[tibble]{tibble}} with information
#'  about data points.
#'
#' @seealso \code{\link[mgcv]{gam}}
#' @importFrom broom augment
#' @importFrom dplyr mutate
#' @importFrom dplyr bind_rows
#' @importFrom mgcv predict.gam
#' @export
#' @examples
#' mean_fit <- NEON_PRIN_5min_cleaned %>%
#'   dplyr::filter(site == "upstream") %>%
#'   dplyr::select(turbidity, level, conductance, temperature) %>%
#'   conditional_mean(turbidity ~ level + conductance + temperature,
#'     knots_mean = c(8, 8, 8)
#'   )
#'
#' data_inf <- mean_fit %>% augment()
augment <- function(x) {
  pred.orig <- predict(x, type = "terms")
  partial.resids <- tibble::as_tibble(pred.orig + residuals(x))

  data <- broom::augment(x) %>%
    dplyr::mutate(
      .cond_EX = as.numeric(mgcv::predict.gam(x)),
      partial.resids
    )

  return(data)
}
