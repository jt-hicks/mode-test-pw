#' Run named odin model

run_model <- function(model = "odin_model",
                           het_brackets = 5,
                           age = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3.5,5,7.5,10,15,20,30,40,50,60,70,80),
                           init_EIR = 10,
                           init_ft = 0.4,
                           country = NULL,
                           admin2 = NULL,
                           time = 100,
                           roster_file = NULL,
                           ...){

  ## create model param list using necessary variables
  mpl <- model_param_list_create(...)

  # generate initial state variables from equilibrium solution
  state <- equilibrium_init_create(age_vector=age, EIR=init_EIR,ft=init_ft,
                                   model_param_list = mpl, het_brackets=het_brackets,
                                   rA_preg = rA_preg
                                   #,
                                   #country = country,
                                   #admin_unit = admin2
                                   )

  # create odin generator
  generator <- switch(model,
    "odin_model" = odin_model,
    "odin_model_emanators" = odin_model_emanators,
    "odin_model_hrp2" = odin_model_hrp2,
    "odin_model_IVM_SMChet" = odin_model_IVM_SMChet,
    "odin_model_TBV" = odin_model_TBV,
    stop(sprintf("Unknown model '%s'", model)))

  # There are many parameters used that should not be passed through
  # to the model.
  state_use <- state[names(state) %in% coef(generator)$name]
  # create model with initial values
  mod <- generator(user = state_use, use_dde = TRUE)
  tt <- seq(0, time, 1)
  # run model
  n_particles <- 10
  start_time <- Sys.time()
  for (j in seq_len(n_particles)) {
    mod_run <- mod$run(tt, verbose=FALSE,step_size_max=9)
  }
  end_time <- Sys.time()
  message(end_time - start_time)

  # shape output
  out <- mod$transform_variables(mod_run)

  # return mod
  return(out)
}
