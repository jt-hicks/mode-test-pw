library("odin.dust")
library("odin")
library("patchwork")
source("MiP-given/model_parameters.R")
source("MiP-given/equilibrium-init-create-stripped.R")

init_age <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)

init_EIR <- 100
prop_treated <- 0.4
rA_preg <- 0.00512821
rU_preg <- 0.00906627
het_brackets <- 5

################## generate the data ######################
# generate random walk of EIR (recursive fn)
genRandWalk <- function(x,vol,randWalk) {
  if (x == 0)    return (randWalk)
  else return(genRandWalk(x-1,vol,c(randWalk,exp(log(randWalk[length(randWalk)])+rnorm(1)*vol))))
}
EIR_times<-seq(0,1800,by=30)

### just a random walk on logscale
EIR_volatility<-0.05
EIR_vals=genRandWalk(length(EIR_times)-1,EIR_volatility,init_EIR)

##set up the simulation for the simualted data 
time<- 5*365
out_step=30

mpl <- model_param_list_create(init_EIR = init_EIR,
                               init_ft = prop_treated,
                               EIR_times=EIR_times,
                               EIR_vals=EIR_vals
                               )

pars <- equilibrium_init_create_stripped(age_vector = init_age,
                                EIR = init_EIR,
                                ft = prop_treated,
                                model_param_list = mpl,
                                het_brackets = het_brackets)

##The malaria model but only on human side (i.e. no mosquitoes to worry about)
generator <- odin("original_malaria/odin_model_stripped_matched.R")
state_use <- pars[names(pars) %in% coef(generator)$name]

# create model with initial values
mod <- generator(user = state_use, use_dde = TRUE)
tt <- seq(0, time, out_step)

# run the simulation to base the data
start.time <- Sys.time()
mod_run <- mod$run(tt)
print(Sys.time()-start.time)

# shape output
out <- mod$transform_variables(mod_run)

# plot data and generate data
plot(out$t,out$prev,col="white")
lines(out$t,out$prev,col="blue",lwd=4)
tested<-round(rnorm(length(out$prev),100,30))
positive<-rbinom(length(out$prev),tested,out$prev)
data_raw<-data.frame(t=out$t+30,tested=tested,positive=positive)
#######################################

######## run particle filter with same model with same log(EIR) random walk but within odin.dust
data <- mcstate::particle_filter_data(data_raw, time = "t", rate = NULL, initial_time = 0)


compare <- function(state, observed, pars = NULL) {
  dbinom(x = observed$positive,
         size = observed$tested,
         prob = state[1, ],
         log = TRUE)
}

index <- function(info) {
  list(run = c(inc = info$index$prev),
       state = c(prev = info$index$prev))
}
stochastic_schedule <- seq(from = 60, by = 30, to = 1830)
#### NB the volatility and initial EIR is hard-coded in the odinmodelmatchedstoch bw lines 230 and 234###
model <- odin.dust::odin_dust("original_malaria/odinmodelmatchedstoch.R")
n_particles <- 100
set.seed(1)

### single with no parallelisation
p_single <- mcstate::particle_filter$new(data, model, n_particles, compare,
                                  index = index, seed = 1L,
                                  stochastic_schedule = stochastic_schedule
                                  )

### takes about 15-20 seconds on my desktop 
start.time <- Sys.time()
lik_single <- p_single$run(pars, save_history = TRUE)
print(Sys.time()-start.time)
history_single <- p_single$history()

###but looks absolutely lauurrrvely
matplot(data_raw$t, t(history_single[1, , -1]), type = "l",
        xlab = "Time", ylab = "State",
        col = "#ff000022", lty = 1, ylim = range(history_single))
lines(out$t,out$prev,col="blue",lwd=4)

## reset seed
set.seed(1)
## run over 4 threads
p_para <- mcstate::particle_filter$new(data, model, n_particles, compare,
                                         index = index, seed = 1L,
                                         stochastic_schedule = stochastic_schedule,
                                       n_threads=4
)
## now takes 3-5 secs
start.time <- Sys.time()
lik_para<- p_para$run(pars, save_history = TRUE)
print(Sys.time()-start.time)

history_para <- p_para$history()
## but looks a bit less perfect :(
matplot(data_raw$t, t(history_para[1, , -1]), type = "l",
        xlab = "Time", ylab = "State",
        col = "#ff000022", lty = 1, ylim = range(history_para))
lines(out$t,out$prev,col="blue",lwd=4)