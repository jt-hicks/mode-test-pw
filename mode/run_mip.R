#Call functions housed in other files for simplicity
#Function that calculates the initial equilibrium state
source("mode/equilibrium-init-create.R")
#Function that defines model parameters
source("mode/model_parameters.R")

# create a vector of age categories
init_age <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)

# provide a value of the annual EIR for this model run
init_EIR <- 10

# provide the length of time (in days) that you want to run the model for
time_period <- 100 #365 * 5

# provide a value for the proportion of cases that are treated (referred to as ft in the paper)
prop_treated <- 0.4

#Create beta (mosquito emergence parameter) sequence to pass to model
betaa_vals <- seq(from = 0.65, to = 3.5, by = 0.15)
betaa_times <- seq(0, time_period, by = ((time_period) / (20 - 1)))

# pregnancy category initial values
rA_preg <- 0.00512821 #same as main compartments
#rA_preg <- 0.1
rU_preg <- 0.00906627 #same as main compartments

het_brackets <- 5
time <- 100
## create model param list using necessary variables
mpl <- model_param_list_create(age = init_age,
                               init_EIR = init_EIR,
                               init_ft = prop_treated,
                               time = time_period,
                               betaa_vals = betaa_vals,
                               betaa_times = betaa_times,
                               rA_preg = rA_preg, rU_preg = rU_preg,
                               comm_age_min = 0,
                               comm_age_max = 59 / 12,
                               anc_age_min = 15,
                               anc_age_max = 20,
                               lag_rates = 10,
                               lag_ratesMos = 10)

# generate initial state variables from equilibrium solution
state <- equilibrium_init_create(age_vector = init_age,
                                 EIR = init_EIR,
                                 ft = prop_treated,
                                 model_param_list = mpl,
                                 het_brackets = het_brackets,
                                 rA_preg = rA_preg)

# create odin generator
res <- pkgload::load_all("mip")
generator <- res$env$mipodinmodel

# create model with initial values
mod <- generator$new(state, time = 0, n_particles = 1)
# run model
res <- seq_along(betaa_times)
info <- mod$info()
mod$set_index(c(prev = info$index$prev))

for (t in seq_along(betaa_times)){
  res[[t]] <- mod$run(betaa_times[[t]])["prev", ]
  mod$update_state(state = betaa_vals[[t]], index = 1)
}

#plot(betaa_times, res)
