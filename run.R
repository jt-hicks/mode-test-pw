case_compare <- function(state, observed) {
  proportion_modelled <- state[1, ] #grab modelled infectious proportion from particle trajectory
  sample_observed <- observed$tested #grab observed sample size
  positive_observed <- observed$positive #grab number of positive tests
  #Get the log density of the difference between the observed and modelled prevalence
  #based on the poisson distribution
  dbinom(x = positive_observed, size = sample_observed, prob = proportion_modelled, log = TRUE)
}

scale_log_weights <- function(log_weights) {
  log_weights[is.nan(log_weights)] <- -Inf
  max_log_weights <- max(log_weights)
  if (!is.finite(max_log_weights)) {
    ## if all log_weights at a time-step are -Inf, this should
    ## terminate the particle filter and output the marginal
    ## likelihood estimate as -Inf
    average <- -Inf
    weights <- rep(NaN, length(log_weights))
  } else {
    ## calculation of weights, there is some rescaling here to avoid
    ## issues where exp(log_weights) might give computationally zero
    ## values
    weights <- exp(log_weights - max_log_weights)
    average <- log(mean(weights)) + max_log_weights
  }
  list(weights = weights, average = average)
}

observed <- read.csv('casedata_monthly.csv')
n_particles <- 1
beta_volatility <- 0.5
data <- unname(split(observed, seq_len(nrow(observed))))
log_likelihood <- 0

gen <- mode::mode("malaria.cpp")
mod <- gen$new(list(), 0, n_particles)
save_history_index <- 1:5
history_value <- array(NA_real_, c(length(save_history_index), n_particles, length(data)))
history_order <- array(NA_integer_, c(n_particles, length(data)))
mod$set_index(2)
d <- data[[1]]
for (i in seq_along(data)) {
  d <- data[[i]]
  y <- mod$run(d$t)
  weights <- scale_log_weights(case_compare(y, d))
  log_likelihood <- log_likelihood + weights$average
  kappa <- sample.int(n_particles, prob = weights$weights, replace = TRUE)
  history_value[, , i] <- mod$state()[save_history_index, ]
  history_order[, i] <- kappa
  mod$reorder(kappa)
  y <- mod$state()
  y[, 5] <- y[, 5] * exp(rnorm(n_particles) * beta_volatility)
  mod$update_state(state = y, reset_step_size = FALSE)
}
