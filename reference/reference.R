gen <- odin::odin("model.R")

case_compare <- function(state, observed, pars) {
  proportion_modelled <- state[["Ih"]]
  dbinom(x = observed$positive,
         size = observed$tested,
         prob = proportion_modelled,
         log = TRUE)
}

observed <- read.csv("../casedata_monthly.csv")
## Very stripped down version of the mcstate::particle_filter_data
## logic
data <- cbind(t_start = observed$t[-nrow(observed)],
              t_end = observed$t[-1],
              observed[-1, -1])

n_particles <- 100

pars <- list(init_Ih = 0.8,
             init_Sv = 100,
             init_Iv = 1,
             nrates = 15,
             init_beta = -log(0.9)) # must match mu

beta_volatility <- 0.5

## Annoyingly there is really no good way of getting a named input
## vector out of odin (see mrc-3156), so there's a bit of a fight here
## to pull it off.
mod <- gen$new(user = pars)
idx <- seq_along(mod$initial(0)) + 1
y0 <- mod$run(c(0, 1))[1, ]
state <- matrix(y0, length(y0), n_particles,
                dimnames = list(names(y0), NULL))

set.seed(1)
log_likelihood <- 0

for (i in seq_len(nrow(data))) {
  d <- data[i, ]

  log_weight <- numeric(n_particles)
  for (j in seq_len(n_particles)) {
    state[, j] <- mod$run(c(d$t_start, d$t_end), state[idx, j])[2, ]
    log_weight[[j]] <- case_compare(state[, j], d, pars)
  }

  weight <- exp(log_weight)
  log_likelihood <- log_likelihood + log(mean(weight))

  kappa <- sample.int(n_particles, prob = weight, replace = TRUE)
  state <- state[, kappa]
  state["beta", ] <- state["beta", ] * exp(rnorm(n_particles) * beta_volatility)
}
