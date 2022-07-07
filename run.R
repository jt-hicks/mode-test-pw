source("MiP-given/model_parameters.R")
source("MiP-given/equilibrium-init-create.R")

data_raw <- read.csv("casedata_monthly.csv",
                     stringsAsFactors = FALSE)
data <- mcstate::particle_filter_data(data_raw, time = "t", rate = NULL, initial_time = 0)

compare <- function(state, observed, pars = NULL) {
  dbinom(x = observed$positive,
         size = observed$tested,
         prob = state[1, ],
         log = TRUE)
}

index <- function(info) {
  list(run = c(inc = info$index$inc),
       state = c(prev = info$index$prev))
}

stochastic_schedule <- seq(from = 30, by = 30, to = 1800)

model <- odin.dust::odin_dust("mipodinmodel.R")
n_particles <- 100
p <- mcstate::particle_filter$new(data, model, n_particles, compare,
                                  index = index, seed = 1L,
                                  stochastic_schedule = stochastic_schedule)

init_age <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)

init_EIR <- 10
prop_treated <- 0.4
rA_preg <- 0.00512821
rU_preg <- 0.00906627
het_brackets <- 5

mpl <- model_param_list_create(age = init_age,
                               init_EIR = init_EIR,
                               init_ft = prop_treated,
                               rA_preg = rA_preg,
                               rU_preg = rU_preg,
                               comm_age_min = 0,
                               comm_age_max = 59 / 12,
                               anc_age_min = 15,
                               anc_age_max = 20,
                               lag_rates = 10,
                               lag_ratesMos = 10)

# generate initial state variables from equilibrium solution
pars <- equilibrium_init_create(age_vector = init_age,
                                 EIR = init_EIR,
                                 ft = prop_treated,
                                 model_param_list = mpl,
                                 het_brackets = het_brackets,
                                 rA_preg = rA_preg)

lik <- p$run(pars, save_history = TRUE)
history <- p$history()
matplot(data_raw$t, t(history[1, , -1]), type = "l",
        xlab = "Time", ylab = "State",
        col = "#ff000022", lty = 1, ylim = range(history))
