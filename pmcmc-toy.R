library("odin.dust")
library("odin")
library("patchwork")
library('mcstate')
library(RColorBrewer)
source('data_gen_toy.R')
#Generate some data
data_raw_toy <- data_gen_toy(0.1)
data <- mcstate::particle_filter_data(data_raw_toy, time = "t", rate = NULL, initial_time = 0)

index <- function(info) {
  list(run = c(Ih = info$index$Ih),
       state = c(Host_prev = info$index$Host_prev,
                 Vector_prev = info$index$Vector_prev))
}

compare <- function(state, observed, pars = NULL) {
  dbinom(x = observed$positive,
         size = observed$tested,
         prob = state[1, ],
         log = TRUE)
}

stochastic_schedule <- seq(from = 60, by = 30, to = 1830)

model <- odin.dust::odin_dust("toyodinmodel.R")
n_particles <- 100
p <- mcstate::particle_filter$new(data, model, n_particles, compare,
                                  index = index, seed = 1L,
                                  stochastic_schedule = stochastic_schedule,
                                  n_threads = 4)

init_Ih <- mcstate::pmcmc_parameter("init_Ih", 0.8, min = 0, max = 1)
init_Sv <- mcstate::pmcmc_parameter("init_Sv", 100, min = 0)
init_Iv <- mcstate::pmcmc_parameter("init_Iv", 1, min = 1)
nrates <- mcstate::pmcmc_parameter("nrates", 15, min = 1)
beta_volatility <- mcstate::pmcmc_parameter("beta_volatility", 0.3, min = 0)

lik <- p$run(pars, save_history = TRUE)
history <- p$history()
matplot(data_raw$t, t(history[1, , -1]), type = "l",
        xlab = "Time", ylab = "State",
        col = "#ff000022", lty = 1)
lines(data_raw$t,data_raw$positive/data_raw$tested,col="blue",lwd=4)

proposal_matrix <- diag(0.1, 5)
mcmc_pars <- mcstate::pmcmc_parameters$new(list(init_Ih = init_Ih, 
                                                init_Sv = init_Sv,
                                                init_Iv = init_Iv,
                                                nrates = nrates,
                                                beta_volatility = beta_volatility),
                                           proposal_matrix)
n_steps<- 500
n_burnin <- 200


control <- mcstate::pmcmc_control(
  n_steps,
  save_state = TRUE,
  save_trajectories = TRUE,
  progress = TRUE,
  n_chains = 4,
  n_workers = 2,
  n_threads_total = 16)
pars<- mcmc_pars$fix(list(init_Sv = 100,
                          init_Iv = 1,
                          nrates = 15))
pars$initial()
pmcmc_run <- mcstate::pmcmc(pars, p, control = control)

history <- pmcmc_run$trajectories$state

matplot(data_raw$t, t(history[1, , -1]), type = "l",
        xlab = "Time", ylab = "State",
        col = "#ff000022", lty = 1)
lines(data_raw$t,data_raw$positive/data_raw$tested,col="blue",lwd=4)

out_pars <- as.data.frame(pmcmc_run$pars)
plot(out_pars$beta_volatility)
plot(out_pars$beta_volatility)
