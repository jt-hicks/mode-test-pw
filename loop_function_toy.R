pf_loop <- function(data_file,n_loop,n_particles=100,volatility=0.5,freq=1){
  data_raw <- read.csv(data_file,
                       stringsAsFactors = FALSE)
  data <- mcstate::particle_filter_data(data_raw, time = "t", rate = NULL, initial_time = 0)
  
  index <- function(info) {
    list(run = c(Ih = info$index$Ih),
         state = c(Host_prev = info$index$Host_prev))
  }
  
  compare <- function(state, observed, pars = NULL) {
    dbinom(x = observed$positive,
           size = observed$tested,
           prob = state[1, ],
           log = TRUE)
  }
  
  lik_list <- function(i,p,pars){
    p$run(pars)
  }
  stochastic_schedule <- seq(from = 60, by = 30/freq, to = 1830)
  
  model <- odin.dust::odin_dust("toyodinmodel.R")

  pars <- list(init_Ih = 0.8,
               init_Sv = 100,
               init_Iv = 1,
               nrates = 15)
  n_particles <- c(10,50,100,200,500)
  lik_hist <- matrix(nrow=n_loop,ncol=length(n_particles))
  #parallel::clusterExport(cl=NULL,varlist = c('mode_toyodinmodel_alloc'))
  start.time <- Sys.time()
  for (k in 1:length(n_particles)){
    p <- mcstate::particle_filter$new(data, model, n_particles[k], compare,
                                      index = index, seed = 1L,
                                      stochastic_schedule = stochastic_schedule,
                                      n_threads=4)
    lik_hist[,k] <- t(unlist(lapply(1:n_loop,FUN = lik_list, p=p, pars=pars)))
    # lik_hist[,k] <- t(unlist(lapply(1:n_loop,lik_list, p=p, pars=pars)))
  }
  print(Sys.time()-start.time)
  df <- as.data.frame(lik_hist)
  colnames(df) <- n_particles
  return(df)
}
