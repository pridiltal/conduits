#' Estimating time delay between two sensors in a river system
#'
#' This function estimates the time that takes water to flow from an upstream location to a downstream
#' location conditional on the observed water-quality variables from the upstream sensor. That time lag is
#' defined as the lag that gives maximum cross-correlation conditional on upstream water-quality variables.
#'
#' @param x 	Model object of class "conditional_ccf" returned from
#'  \code{\link[conduits]{conditional_ccf}}
#' @return A \code{\link[tibble]{tibble}} with estimated time lag "dt"
#'  and corresponding maximum cross-correlation
#' @importFrom dplyr mutate
#' @importFrom stats predict.glm pnorm
#' @importFrom purrr map_dfc
#'
#' @author Puwasala Gamakumara & Priyanga Dilini Talagala
#'
#' @examples
#'
#'old_ts <- NEON_PRIN_5min_cleaned %>%
#'  dplyr::select(
#'    Timestamp, site, turbidity, level,
#'    conductance, temperature
#'  ) %>%
#' tidyr::pivot_wider(
#'    names_from = site,
#'    values_from = turbidity:temperature
#'  )
#'
#' fit_mean_y <- old_ts %>%
#'   conditional_mean(turbidity_downstream ~
#'                      s(level_upstream, k = 8) +
#'                      s(conductance_upstream, k = 8) +
#'                      s(temperature_upstream, k = 8))
#'
#' fit_var_y <- old_ts %>%
#'   conditional_var(
#'     turbidity_downstream ~
#'       s(level_upstream, k = 7) +
#'       s(conductance_upstream, k = 7) +
#'       s(temperature_upstream, k = 7),
#'     family = "Gamma",
#'     fit_mean_y
#'   )
#'
#' fit_mean_x <- old_ts %>%
#'   conditional_mean(turbidity_upstream ~
#'                      s(level_upstream, k = 8) +
#'                      s(conductance_upstream, k = 8) +
#'                      s(temperature_upstream, k = 8))
#'
#' fit_var_x <- old_ts %>%
#'   conditional_var(
#'     turbidity_upstream ~
#'       s(level_upstream, k = 7) +
#'       s(conductance_upstream, k = 7) +
#'       s(temperature_upstream, k = 7),
#'     family = "Gamma",
#'     fit_mean_x
#'  )
#'
#' fit_c_ccf <- old_ts %>%
#'    tidyr::drop_na() %>%
#'    conditional_ccf(
#'      I(turbidity_upstream*turbidity_downstream) ~ splines::ns(
#'      level_upstream, df = 5) +
#'      splines::ns(conductance_upstream, df = 5),
#'      lag_max = 10,
#'      fit_mean_x, fit_var_x, fit_mean_y, fit_var_y,
#'      df_correlation = c(5,5))
#'
#' new_data <- fit_c_ccf %>% estimate_dt()
#'
#' @export
#'
estimate_dt <- function(x){

  data_inf <- x %>% augment()
  k <- x$lag_max
  cnames <- paste("c", 1:k, sep = "")
  data_sub <- data_inf %>%
    dplyr::mutate(dt =  max.col(.[cnames], ties.method="first"))

  data_sub$max_ccf <- apply(data_sub[cnames],1,max)


  # p value calculation
  stdz_data_predict_ccf <- matrix(0, ncol = length(k), nrow = nrow(data_sub))

  predict_ccf_gam_stdz <- function(t)
  {
    cond_ccf <- stats::predict.glm(x[[t]], newdata = x$data, type = "response",
                                   se.fit = TRUE)
    cond_ccf_std <- cond_ccf$fit/ cond_ccf$se.fit
    return(cond_ccf_std)
  }

  c_ccf_est_std <- purrr::map_dfc(1:k, predict_ccf_gam_stdz)
  row_max <- apply(c_ccf_est_std, 1, max)
  data_sub$pval <- 1 - (stats::pnorm(row_max))^k

  return(data_sub)

}
