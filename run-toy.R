data_raw <- read.csv("casedata_monthly.csv",
                     stringsAsFactors = FALSE)
data <- mcstate::particle_filter_data(data_raw, time = "t", rate = NULL, initial_time = 0)

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

stochastic_schedule <- seq(from = 30, by = 30, to = 1800)

model <- odin.dust::odin_dust("toyodinmodel.R")
n_particles <- 100
p <- mcstate::particle_filter$new(data, model, n_particles, compare,
                                  index = index, seed = 1L,
                                  stochastic_schedule = stochastic_schedule)

pars <- list(init_Ih = 0.8,
             init_Sv = 100,
             init_Iv = 1,
             nrates = 15)

lik <- p$run(pars, save_history = TRUE)
history <- p$history()
matplot(data_raw$t, t(history[1, , -1]), type = "l",
        xlab = "Time", ylab = "State",
        col = "#ff000022", lty = 1, ylim = c(0, 1))
matlines(data_raw$t, t(history[2, , -1]), col = "#0000ff22", lty = 1)
