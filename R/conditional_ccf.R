#' Computing conditional cross-correlations at given lags
#'
#' This function computes cross correlation between $x_t$ and $y_{t+k}$ at $k = 1,2,...$
#' conditional on a set of time series $z_t$
#'
#' @param data a tibble containing all the time series including
#' ystar*xstar which are uniquely identified by the corresponding
#' Timestamp.
#' @param formula A GAM formula. See \code{\link[mgcv]{formula.gam}}.
#' @param lag_max Maximum lag at which to calculate the conditional ccf
#' @param df_correlation a vector specifying the degrees of freedom to be considered for each numerical
#' predictor when fitting additive models for conditional cross-correlations. Each component of the
#' vector should corresponds to each predictor specified in "z_numeric".
#' @return The function returns a list of objects of class
#' "glm" as described in \code{\link[stats]{glm}}.
#'
#' @details{ Suppose $x_t$ and $y_t$ are conditionally normalised with respect
#' to $z_t$ using \code{conditional_mean} and \code{conditional_var}. Then
#' we can estimate the conditional cross-correlation between $x_t$ and $y_t$ at lag $k$, i.e. $r_k = E(x_ty_{t+k}|z_t)$
#' via generalised additive models (GAM). \code{conditional_ccf} uses natural splines implemented
#' in \code{splines} package to estimate the conditional cross-correlations between two
#' time series given a set of time series predictors. Users first need  to
#' normalise $x_t$ and $y_t$ at lag $k$ using \code{conditional_mean} and \code{conditional_var}}
#'
#' @seealso \code{\link[stats]{glm}}
#' @importFrom stats update glm as.formula
#' @importFrom purrr map
#'
#' @export
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
#' calc_xyk_star <- function(k, old_ts, fit_mean_y, fit_var_y) {
#'   old_ts_lead <- old_ts %>%
#'     dplyr::mutate_at("turbidity_downstream",
#'                      dplyr::lead,
#'                      n = k
#'     ) %>%
#'     normalize(
#'       ., turbidity_downstream,
#'       fit_mean_y,
#'       fit_var_y
#'     ) * old_ts$xstar
#' }
#'
#' k <- 3
#'
#' new_ts <- old_ts %>%
#'  dplyr::mutate(xstar = normalize(
#'     ., turbidity_upstream,
#'    fit_mean_x, fit_var_x
#'   )) %>%
#'   purrr::map_dfc(
#'     1:k, calc_xyk_star, .,
#'     fit_mean_y, fit_var_y
#'  ) %>%
#'   stats::setNames(paste("xystar_t", 1:k, sep = "")) %>%
#'   dplyr::bind_cols(old_ts, .)
#'
#' fit_c_ccf <- new_ts %>%
#'    tidyr::drop_na() %>%
#'    conditional_ccf(
#'      xystar ~ splines::ns(
#'      level_upstream, df = 5) +
#'      splines::ns(conductance_upstream, df = 5),
#'      lag_max = k,
#'      df_correlation = c(5,5))
#'
#'
conditional_ccf <- function(data, formula, lag_max, df_correlation){

  vars <- all.vars(formula)
  y_name <- vars[1]
  xynames <- colnames(data)[grepl("xystar" , names(data ))]
  corrl <- corrlink()

  fit_ccf_gam <- function(k)
  {
    fk<- paste(xynames[k], "~.")
    formula_k <- stats::update(formula,  stats::as.formula(fk))
    ccf_gam_fit_k <- stats::glm(formula = formula_k,
                                   data = data,
                                   family = stats::gaussian(link = corrl),
                                   start = rep(0,(sum(df_correlation)+1)),
                                   control = stats::glm.control(maxit = 400))
    return(ccf_gam_fit_k)
  }

  ccf_gam_fit <- purrr::map(1:lag_max, fit_ccf_gam)
  ccf_gam_fit$data <- data
  class(ccf_gam_fit) <- c("conditional_ccf")
  return(ccf_gam_fit)
}

corrlink <- function() {
  ## link
  linkfun <- function(mu) {log((1+mu)/(1-mu))}
  ## inverse link
  linkinv <- function(eta) {(exp(eta) - 1)/(exp(eta) + 1)}
  ## derivative of invlink wrt eta
  mu.eta <- function(eta) { 2*exp(eta)/(exp(eta) + 1)^2 }
  valideta <- function(eta) TRUE
  link <- "corrlink"
  structure(list(linkfun = linkfun,
                 linkinv = linkinv,
                 mu.eta = mu.eta,
                 valideta = valideta,
                 name = link),
            class = "link-glm")
}
