remotes::install_github("mrc-ide/odin.dust", upgrade = TRUE)
remotes::install_github("mrc-ide/mcstate", upgrade = TRUE, force = TRUE)

library("odin.dust")
library("odin")
library("patchwork")
library('mcstate')
library(RColorBrewer)
source("MiP-given/model_parameters.R")
source("MiP-given/equilibrium-init-create-stripped.R")
library("coda")

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
EIR_volatility<-0.4
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
saveRDS(data_raw,'data_raw_120822.RDS')
data_raw <- readRDS('data_raw_120822.RDS')
########################################
####### Function to plot particle trajectories
plot_particle_filter <- function(history, true_history, times, obs_end = NULL) {
  if (is.null(obs_end)) {
    obs_end <- max(times)
  }
  true_history$prev <- true_history$positive/true_history$tested
  par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
  matplot(times, t(history[1, , -1]), type = "l",
          xlab = "Time", ylab = "Prevalence",
          col = "#A6CEE3", lty = 1, ylim = range(history))
  matpoints(times, true_history$prev, pch = 19,
            col = "#1F78B4")
}


#######################################

######## run pMCMC with same model with same log(EIR) random walk but within odin.dust

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
mpl_pf <- model_param_list_create(EIR_SD=0.4)

state <- equilibrium_init_create_stripped(age_vector = init_age,
                                            EIR = init_EIR,
                                            ft = prop_treated,
                                            model_param_list = mpl_pf,
                                            het_brackets = het_brackets)
#### NB the volatility and initial EIR is hard-coded in the odinmodelmatchedstoch bw lines 230 and 234###
model <- odin.dust::odin_dust("original_malaria/odinmodelmatchedstoch.R")


n_particles <- 100
set.seed(1)

### Set particle filter
pf <- mcstate::particle_filter$new(data, model, n_particles, compare,
                                         index = index, seed = 1L,
                                         stochastic_schedule = stochastic_schedule,
                                         ode_control = mode::mode_control(max_steps = 1e7),
                                         n_threads = 4)
start.time <- Sys.time()
lik <- pf$run(state, save_history = TRUE)
print(Sys.time()-start.time)
history <- pf$history()
plot_particle_filter(history,true_history=out,times=out$t)

### Set pmcmc control
n_steps <- 500
n_burnin <- 200
control <- mcstate::pmcmc_control(
  n_steps,
  save_state = TRUE,
  save_trajectories = TRUE,
  progress = TRUE,
  n_threads_total = 16)

### Set pmcmc parameters
##Need to convert equilibrium into pmcmc parameters?
pmcmc_par <- function(name){
  print(name)
  mcstate::pmcmc_parameter(name,initial = 1)

}
EIR_SD <- mcstate::pmcmc_parameter("EIR_SD", 0.3, min = 0)
equil_list <- lapply(names(state),pmcmc_par)
names(equil_list) <- names(state)
proposal_matrix <- diag(0.1, length(equil_list))
list <- append(list(EIR_SD = EIR_SD), equil_list[-length(equil_list)])
mcmc_pars <- mcstate::pmcmc_parameters$new(list,
                                           proposal_matrix)
mcmc_pars <- mcmc_pars$fix(state[-length(state)])
pmcmc_parameters$
### Run pMCMC
start.time <- Sys.time()
pmcmc_run <- mcstate::pmcmc(mcmc_pars, pf, control = control)
print(Sys.time()-start.time)
history <- pmcmc_run$trajectories$state
saveRDS(pmcmc_run,'C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PMCMC/Development/pmcmc/pmcmc_run_120822.rds')
plot_particle_filter(history,true_history=data_raw,times=data_raw$t)
pars_post <- as.data.frame(pmcmc_run$pars)
plot(pars_post$EIR_SD)
ll <- as.data.frame(pmcmc_run$probabilities)
plot(ll$log_likelihood)
mcmc <- coda::as.mcmc(cbind(
  pmcmc_run$probabilities, pmcmc_run$pars))
1 - coda::rejectionRate(mcmc)
coda::effectiveSize(mcmc)
cov(pmcmc_run$pars)
summary(mcmc)
windows(60,50)
plot(mcmc)
