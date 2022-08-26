pf_loop <- function(volatility=0.5,data_file,n_loop,n_particles=200,freq=1,
                    init_EIR=100,atol=NULL,rtol=NULL,max_steps=NULL){
  start.time <- Sys.time()

  init_age <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)
  prop_treated <- 0.4
  rA_preg <- 0.00512821
  rU_preg <- 0.00906627
  het_brackets <- 5

  data <- mcstate::particle_filter_data(as.data.frame(data_file), time = "t", rate = NULL, initial_time = 0)
  
  index <- function(info) {
    list(run = c(inc = info$index$prev),
         state = c(prev = info$index$prev))
  }
  
  compare <- function(state, observed, pars = NULL) {
    dbinom(x = observed$positive,
           size = observed$tested,
           prob = state[1, ],
           log = TRUE)
  }
  
  lik_list <- function(i,p,pars){
    start.time <- Sys.time()
    run <- tryCatch(p$run(pars), 
                    error=function(e) conditionMessage(e))
    time <- difftime(Sys.time(), start.time, units = "secs")
    temp <- c(likelihood=run,time=time)
    return(temp)
  }
  stochastic_schedule <- seq(from = 60, by = 30/freq, to = 1830)
  
  model <- odin.dust::odin_dust("original_malaria/odinmodelmatchedstoch.R")
  
  mpl <- model_param_list_create(EIR_SD = volatility)
  
  pars <- equilibrium_init_create_stripped(age_vector = init_age,
                                           EIR = init_EIR,
                                           ft = prop_treated,
                                           model_param_list = mpl,
                                           het_brackets = het_brackets)
  lik_hist <- matrix(nrow=n_loop,ncol=length(n_particles))
  time_hist <- matrix(nrow=n_loop,ncol=length(n_particles))
  print(Sys.time()-start.time)
  start.time <- Sys.time()

  for (k in 1:length(n_particles)){
    p <- mcstate::particle_filter$new(data, model, n_particles[k], compare,
                                      index = index, #seed = 1L,
                                      stochastic_schedule = stochastic_schedule,
                                      ode_control = mode::mode_control(max_steps=max_steps, atol=atol, rtol=rtol),
                                      n_threads=4)
    # lik_hist[,k] <- t(unlist(parallel::parLapply(NULL,1:n_loop,fun = lik_list, p=p, pars=pars)))
    temp <- data.frame(likelihood=numeric(),time=numeric())
    temp <- lapply(1:n_loop,lik_list, p=p, pars=pars)
    lik_hist[,k] <- t(unlist(sapply(temp,'[',1)))
    time_hist[,k] <- t(unlist(sapply(temp,'[',2)))
  }
  print(Sys.time()-start.time)
  df1 <- as.data.frame(lik_hist)
  df1$var <- 'likelihood'
  df2 <- as.data.frame(time_hist)
  df2$var <- 'time'
  df <- rbind(df1,df2)
 colnames(df) <- c(n_particles,'variable')
  return(df)
}
