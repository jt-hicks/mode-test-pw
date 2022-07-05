case_compare <- function(state, observed) {
  proportion_modelled <- state[1, ]
  dbinom(x = observed$positive,
         size = observed$tested,
         prob = proportion_modelled,
         log = TRUE)
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

observed <- read.csv("casedata_monthly.csv")
n_particles <- 100
beta_volatility <- 0.5
data <- unname(split(observed, seq_len(nrow(observed))))[-1]

gen <- mode::mode("malaria.cpp")

history_index <- setNames(1:5, c("Sh", "Ih", "Sv", "Iv", "beta"))
beta_index <- 5L
n_state <- length(history_index)
n_data <- length(data)
history_value <- array(NA_real_, c(n_state, n_particles, n_data))
history_order <- array(NA_integer_, c(n_particles, n_data))

mod <- gen$new(list(), 0, n_particles)
mod$set_index(2)
mod$set_stochastic_schedule(observed$t[-1])

mod$update_state(time = 0)
log_likelihood <- 0
log_likelihood_step <- numeric(n_data)
for (i in seq_len(n_data)) {
  d <- data[[i]]

  y <- mod$run(d$t)

  weights <- scale_log_weights(case_compare(y, d))
  log_likelihood_step[[i]] <- weights$average
  log_likelihood <- log_likelihood + weights$average

  kappa <- sample.int(n_particles, prob = weights$weights, replace = TRUE)

  history_value[, , i] <- mod$state(history_index)
  history_order[, i] <- kappa

  mod$reorder(kappa)
}

## This will assemble the full history, but because we sample down to
## a single point at the end it's not very interesting:
history <- mcstate:::history_single(history_value, history_order,
                                    history_index, index_particle = NULL)
matplot(observed$t[-1], t(history["Ih", , ]), type = "l",
        lty = 1, lwd = 0.5, col = "#00000055")
matplot(observed$t[-1], t(history["beta", , ]), type = "l",
        lty = 1, lwd = 0.5, col = "#00000055")
