res <- pkgload::load_all("malaria")
gen <- res$env$odinmodel

case_compare <- function(proportion_modelled, observed) {
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
## Very stripped down version of the mcstate::particle_filter_data
## logic
data <- cbind(t_start = observed$t[-nrow(observed)],
              t_end = observed$t[-1],
              observed[-1, -1])

n_particles <- 100

pars <- list(init_Ih = 0.8,
             init_Sv = 100,
             init_Iv = 1,
             nrates = 15)

set.seed(1)

log_likelihood <- rep(0, 50)

for (s in seq_along(log_likelihood)) {
  mod <- gen$new(pars, 0, n_particles, seed = 1)
  mod$set_index(c(Ih = 2L))
  mod$set_stochastic_schedule(data$t_end)

  for (i in seq_len(nrow(data))) {
    d <- data[i,]

    log_weight <- numeric(n_particles)
    state <- mod$run(d$t_end)

    for (j in seq_len(n_particles)) {
      log_weight[[j]] <- case_compare(state["Ih", j], d)
    }

    scaled_weight <- scale_log_weights(log_weight)
    log_likelihood[[s]] <- log_likelihood[[s]] + scaled_weight$average

    kappa <- sample.int(n_particles, prob = scaled_weight$weights, replace = TRUE)
    mod$reorder(kappa)
  }

}

mean(log_likelihood)
var(log_likelihood) ^ 0.5