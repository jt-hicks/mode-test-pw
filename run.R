data_raw <- read.csv("casedata_monthly.csv",
                     stringsAsFactors = FALSE)
data <- mcstate::particle_filter_data(data_raw, time = "t", rate = NULL, initial_time = 0)

compare <- function(state, observed, pars = NULL) {
  dbinom(x = observed$positive,
         size = observed$tested,
         prob = state[1,],
         log = TRUE)
}

index <- function(info) {
  list(run = c(Ih = info$index$Ih),
       state = c(Ih = info$index$Ih,
                 Sh = info$index$Sh))
}

stochastic_schedule <- seq(from = 30, by = 30, to = 1800)

model <- odin.dust::odin_dust("mipodinmodel.R")
n_particles <- 100
p <- mcstate::particle_filter$new(data, model, n_particles, compare,
                                  index = index, seed = 1L,
                                  stochastic_schedule = stochastic_schedule)
pars <- list(init_Ih = 0.8,
             init_Sv = 100,
             init_Iv = 1,
             nrates = 15)

lik <- p$run(pars)
