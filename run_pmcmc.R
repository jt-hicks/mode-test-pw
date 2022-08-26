
run_pmcmc <- function(data_raw,
                      n_particles=200,
                      EIR_vol,
                      # proposal_dist,
                      init_age = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80),
                      init_EIR = 100,
                      prop_treated = 0.4,
                      rA_preg = 0.00512821,
                      rU_preg = 0.00906627,
                      het_brackets = 5,
                      max_steps = 1e7,
                      atol = 1e-3,
                      rtol = 1e-6){
  ######## run pMCMC with same model with same log(EIR) random walk but within odin.dust
  
  data <- mcstate::particle_filter_data(data_raw, time = "t", rate = NULL, initial_time = 0)
  
  compare <- function(state, observed, pars = NULL) {
    dbinom(x = observed$positive,
           size = observed$tested,
           prob = state[1, ],
           log = TRUE)
  }
  
  index <- function(info) {
    list(run = c(prev = info$index$prev),
         state = c(prev = info$index$prev,
                   EIR = info$index$EIR_out,
                   inc = info$index$incunder5))
  }
  
  stochastic_schedule <- seq(from = 60, by = 30, to = 1830)
  
  mpl_pf <- model_param_list_create(EIR_SD=EIR_vol)
  
  state <- equilibrium_init_create_stripped(age_vector = init_age,
                                            EIR = init_EIR,
                                            ft = prop_treated,
                                            model_param_list = mpl_pf,
                                            het_brackets = het_brackets)
  
  #### NB the volatility and initial EIR is hard-coded in the odinmodelmatchedstoch bw lines 230 and 234###
  model <- odin.dust::odin_dust("original_malaria/odinmodelmatchedstoch.R")
  
  #set.seed(1)
  
  ### Set particle filter
  pf <- mcstate::particle_filter$new(data, model, n_particles, compare,
                                     index = index, seed = 1L,
                                     stochastic_schedule = stochastic_schedule,
                                     ode_control = mode::mode_control(max_steps = max_steps, atol = atol, rtol = rtol),
                                     n_threads = 4)
  
  ### Set pmcmc control
  n_steps <- 500
  control <- mcstate::pmcmc_control(
    n_steps,
    save_state = TRUE,
    save_trajectories = TRUE,
    progress = TRUE,
    n_chains = 1,
    n_workers = 1,
    n_threads_total = 4)
  
  ### Set pmcmc parameters
  ##Need to convert equilibrium into pmcmc parameters?
  pmcmc_par <- function(name){
    print(name)
    mcstate::pmcmc_parameter(name,initial = 1)
  }
  EIR_SD <- mcstate::pmcmc_parameter("EIR_SD", 0.3, min = 0,max=0.9,
                                     prior = function(p) dexp(p, rate = 5, log = TRUE))
  # init_EIR <- mcstate::pmcmc_parameter("init_EIR", 100, min = 0, max = 1000, 
  #                                      prior = function(p) dgamma(p, shape = 1, scale = 0.01, log = TRUE))
  
  equil_list <- lapply(names(state),pmcmc_par)
  names(equil_list) <- names(state)
  proposal_matrix <- diag(0.005, length(equil_list))
  # proposal_matrix[1,1:2] <- proposal_dist[1,]
  # proposal_matrix[2,1:2] <- proposal_dist[2,]
  list <- append(list(EIR_SD = EIR_SD), equil_list[-length(equil_list)])
  mcmc_pars <- mcstate::pmcmc_parameters$new(list,
                                             proposal_matrix)
  mcmc_pars <- mcmc_pars$fix(state[-length(state)])

  ### Run pMCMC
  start.time <- Sys.time()
  pmcmc_run <- mcstate::pmcmc(mcmc_pars, pf, control = control)
  print(Sys.time()-start.time)
  return(pmcmc_run)
}