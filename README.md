# Test repo for running Malaria in Pregnancy model through a particle filter
This repo contains a bunch of test code from various stages of making this work and most files in sub-directories can be ignored. 
The only relevant files now are those in the root directory and [MiP-given/equilibrium-init-create.R](MiP-given/equilibrium-init-create.R), [MiP-given/model_parameters.R](MiP-given/model_parameters.R), and [MiP-given/MZ_multi_rates.rds](MiP-given/MZ_multi_rates.rds) 
for generating the model parameters.

## Installation
Install the latest versions of [odin.dust](https://github.com/mrc-ide/odin.dust) and 
[mcstate](https://github.com/mrc-ide/mcstate/) from GitHub

```bash
remotes::install_github("mrc-ide/odin.dust", upgrade = FALSE)
remotes::install_github("mrc-ide/mcstate", upgrade = FALSE)
```

## Running
### Toy model example

Run and plot particle trajectories by running the `run-toy.R` script. 
This follows these steps:

Prepare data set using `mcstate::particle_filter_data`.
`rate` must be `NULL` here and `initial_time` will most likely be 0 (the first data value must be > 0, 
so the times in the original data provided have been incremented accordingly)

```r
data_raw <- read.csv("casedata_monthly.csv",
                       stringsAsFactors = FALSE)
data <- mcstate::particle_filter_data(data_raw, time = "t", rate = NULL, initial_time = 0)
```

Define an index function for filtering the run state. The first list item, `run` should contain the 
portion of state needed for the likelihood calculation; the second list item, `state` should contain the 
portion of state for which particle history will be saved (see below)

```r
index <- function(info) {
  list(run = c(Ih = info$index$Ih),
       state = c(Ih = info$index$Ih,
                 Sh = info$index$Sh))
}
```

Define a comparison function
```r
compare <- function(state, observed, pars = NULL) {
    Ih <- state[1, ] # as defined by the above index
    dbinom(x = observed$positive,
           size = observed$tested,
           prob = Ih,
           log = TRUE)
}
```

A schedule for running the stochastic updates 
```r
stochastic_schedule <- seq(from = 60, by = 30, to = 1830)
```

Run the model and get a likelihood

```r
model <- odin.dust::odin_dust("toyodinmodel.R")
n_particles <- 100
p <- mcstate::particle_filter$new(data, model, n_particles, compare,
                       index = index, seed = 1L,
                       stochastic_schedule = stochastic_schedule)
pars <- list(init_Ih = 0.8,
             init_Sv = 100,
             init_Iv = 1,
             nrates = 15)
lik <- p$run(pars)
```

To plot particle trajectories, run the model with `save_history = TRUE`:

```r
lik <- p$run(pars, save_history = TRUE)
history <- p$history()
matplot(data_raw$t, t(history[1, , -1]), type = "l",
        xlab = "Time", ylab = "State",
        col = "#ff000022", lty = 1, ylim = range(history))
matlines(data_raw$t, t(history[2, , -1]), col = "#0000ff22", lty = 1)

```

### Malaria in pregnancy model
Run and plot particle trajectories by running the `run.R` script.
This follows exactly the same steps as above. `betaa_td` in this model isn't actually being updated via a stochastic process,
 but just incremented every time step as per the example you gave us. This can be edited in the model code [mipodinmodel.R](mipodinmodel.R) (see the toy model syntax in [toyodinmodel.R](toyodinmodel.R) as a reference).
The parameters for the model have been generated using the exact same scripts you provided us (now found at [MiP-given/model_parameters.R](MiP-given/model_parameters.R)
and [MiP-given/equilibrium-init-create.R](MiP-given/equilibrium-init-create.R))
