gen_run <- function(EIR_volatility,init_EIR,proposal){
  run_pmcmc(data=data_gen(EIR_volatility,init_EIR = init_EIR),
            # EIR_vol=0.3,
            # init_EIR = 5,
            n_particles = 200,
            max_steps = 5e6,
            atol = 1e-5,
            rtol = 1e-5,
            proposal_matrix = proposal,
            n_steps = 10000)
}
