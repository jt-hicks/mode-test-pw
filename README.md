# Test repo for running Malaria in Pregnancy model through a particle filter
This repo contains a bunch of test code from various stages of making this work. 
All sub-directories can be completely ignored; the final version of the code
is contained in the root directory.

## Installation
Install the development versions of [odin.dust](https://github.com/mrc-ide/odin.dust) and 
[mcstate](https://github.com/mrc-ide/mcstate/) from GitHub

```bash
remotes::install_github("mrc-ide/odin.dust")
remotes::install_github("mrc-ide/mcstate")
```

## Running

Prepare data set using `mcstate::particle_filter_data`.
`rate` must be `NULL` here and `initial_time` will most likely be 0 (the first data value must be > 0, 
so the times in the original data provided have been incremented accordingly)


```r
data_raw <- read.csv("casedata_monthly.csv",
                       stringsAsFactors = FALSE)
data <- mcstate::particle_filter_data(data_raw, time = "t", rate = NULL, initial_time = 0)
```

Define a comparison function
```r
compare <- function(state, observed, pars = NULL) {
    dbinom(x = observed$positive,
           size = observed$tested,
           prob = state[1, ],
           log = TRUE)
}
```

Define an index function for filtering the run state
```r
index <- function(info) {
  list(run = c(Ih = info$index$Ih),
       state = c(Ih = info$index$Ih,
                 Sh = info$index$Sh))
}
```

A schedule for running the stochastic updates 
```r
stochastic_schedule  <- seq(from = 30, by = 30, to = 1800)
```

Run the model and get a likelihood

```r
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
```
