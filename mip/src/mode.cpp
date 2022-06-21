#include <mode/r/mode.hpp>

template <typename real_type, typename container>
__host__ __device__ real_type odin_sum1(const container x, size_t from, size_t to);
template <typename real_type, typename container>
__host__ __device__ real_type odin_sum2(const container x, int from_i, int to_i, int from_j, int to_j, int dim_x_1);
template <typename real_type, typename container>
__host__ __device__ real_type odin_sum3(const container x, int from_i, int to_i, int from_j, int to_j, int from_k, int to_k, int dim_x_1, int dim_x_12);
// [[dust::class(mipodinmodel)]]
// [[dust::param(aD, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(age, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(age_20_factor, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(age_rate, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(age05, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(age20l, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(age20u, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(age59, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(ageend, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(agestart, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(alphaA, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(alphaU, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(av0, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(b0, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(b1, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(cD, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(cT, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(cU, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(d1, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(dB, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(dCA, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(dCM, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(dE, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(delayGam, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(delayMos, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(den, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(dID, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(DY, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(eta, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(fD0, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(foi_age, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(FOI_eq, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(FOIv_eq, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(ft, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(gamma1, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(gammaD, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(het_wt, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(IB0, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(IC0, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(ID0, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(init_A, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(init_A_preg, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(init_D, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(init_D_preg, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(init_Ev, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(init_IB, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(init_ICA, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(init_ICM, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(init_ID, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(init_Iv, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(init_P, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(init_P_preg, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(init_pregs, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(init_S, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(init_S_preg, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(init_Sv, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(init_T, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(init_T_preg, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(init_U, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(init_U_preg, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(kB, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(kC, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(kD, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(lag_rates, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(lag_ratesMos, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(mu0, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(mv0, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(na, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(nh, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(nrates, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(omega, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(p10, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(p2, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(phi0, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(phi1, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(PM, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(rA, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(rA_preg, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(rD, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(rel_foi, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(rP, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(rT, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(rU, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(rU_preg, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(sample_rates, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(sample_transition_rates, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(tau1, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(tau2, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(uB, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(uCA, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(uD, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(wane_rates, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(x_I, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
class mipodinmodel {
public:
  using real_type = double;
  using rng_state_type = dust::random::generator<real_type>;
  using data_type = mode::no_data;
  struct shared_type {
    real_type aD;
    std::vector<real_type> age;
    real_type age_20_factor;
    std::vector<real_type> age_rate;
    int age05;
    int age20l;
    int age20u;
    int age59;
    real_type ageend;
    real_type agestart;
    real_type alphaA;
    real_type alphaU;
    real_type av0;
    real_type b0;
    real_type b1;
    real_type cD;
    real_type cT;
    real_type cU;
    real_type d1;
    real_type dB;
    real_type dCA;
    real_type dCM;
    real_type dE;
    real_type delayGam;
    real_type delayMos;
    std::vector<real_type> den;
    real_type dID;
    int dim_A;
    int dim_A_1;
    int dim_A_2;
    int dim_A_preg;
    int dim_A_preg_1;
    int dim_A_preg_2;
    int dim_age;
    int dim_age_rate;
    int dim_b;
    int dim_b_1;
    int dim_b_2;
    int dim_cA;
    int dim_cA_1;
    int dim_cA_2;
    int dim_clin_inc;
    int dim_clin_inc_1;
    int dim_clin_inc_2;
    int dim_clin_inc_preg;
    int dim_clin_inc_preg_1;
    int dim_clin_inc_preg_2;
    int dim_clin_inc0to5;
    int dim_clin_inc0to5_1;
    int dim_clin_inc0to5_2;
    int dim_D;
    int dim_D_1;
    int dim_D_2;
    int dim_D_preg;
    int dim_D_preg_1;
    int dim_D_preg_2;
    int dim_den;
    int dim_EIR;
    int dim_EIR_1;
    int dim_EIR_2;
    int dim_fd;
    int dim_FOI;
    int dim_FOI_1;
    int dim_FOI_12;
    int dim_FOI_2;
    int dim_FOI_3;
    int dim_foi_age;
    int dim_FOI_eq;
    int dim_FOI_eq_1;
    int dim_FOI_eq_2;
    int dim_FOI_lag;
    int dim_FOI_lag_1;
    int dim_FOI_lag_2;
    int dim_FOIv;
    int dim_FOIvijk;
    int dim_FOIvijk_1;
    int dim_FOIvijk_2;
    int dim_het_wt;
    int dim_IB;
    int dim_IB_1;
    int dim_IB_2;
    int dim_IC;
    int dim_IC_1;
    int dim_IC_2;
    int dim_ICA;
    int dim_ICA_1;
    int dim_ICA_2;
    int dim_ICM;
    int dim_ICM_1;
    int dim_ICM_2;
    int dim_ID;
    int dim_ID_1;
    int dim_ID_2;
    int dim_ince_delay;
    int dim_init_A;
    int dim_init_A_1;
    int dim_init_A_2;
    int dim_init_A_preg;
    int dim_init_A_preg_1;
    int dim_init_A_preg_2;
    int dim_init_D;
    int dim_init_D_1;
    int dim_init_D_2;
    int dim_init_D_preg;
    int dim_init_D_preg_1;
    int dim_init_D_preg_2;
    int dim_init_FOI;
    int dim_init_FOI_1;
    int dim_init_FOI_12;
    int dim_init_FOI_2;
    int dim_init_FOI_3;
    int dim_init_IB;
    int dim_init_IB_1;
    int dim_init_IB_2;
    int dim_init_ICA;
    int dim_init_ICA_1;
    int dim_init_ICA_2;
    int dim_init_ICM;
    int dim_init_ICM_1;
    int dim_init_ICM_2;
    int dim_init_ICM_pre;
    int dim_init_ID;
    int dim_init_ID_1;
    int dim_init_ID_2;
    int dim_init_P;
    int dim_init_P_1;
    int dim_init_P_2;
    int dim_init_P_preg;
    int dim_init_P_preg_1;
    int dim_init_P_preg_2;
    int dim_init_pregs;
    int dim_init_S;
    int dim_init_S_1;
    int dim_init_S_2;
    int dim_init_S_preg;
    int dim_init_S_preg_1;
    int dim_init_S_preg_2;
    int dim_init_T;
    int dim_init_T_1;
    int dim_init_T_2;
    int dim_init_T_preg;
    int dim_init_T_preg_1;
    int dim_init_T_preg_2;
    int dim_init_U;
    int dim_init_U_1;
    int dim_init_U_2;
    int dim_init_U_preg;
    int dim_init_U_preg_1;
    int dim_init_U_preg_2;
    int dim_P;
    int dim_P_1;
    int dim_P_2;
    int dim_p_det;
    int dim_p_det_1;
    int dim_p_det_2;
    int dim_P_preg;
    int dim_P_preg_1;
    int dim_P_preg_2;
    int dim_phi;
    int dim_phi_1;
    int dim_phi_2;
    int dim_pregs;
    int dim_prev_cba;
    int dim_prev_cba_1;
    int dim_prev_cba_2;
    int dim_prev_rdt;
    int dim_prev_rdt_1;
    int dim_prev_rdt_2;
    int dim_prev0to59;
    int dim_prev0to59_1;
    int dim_prev0to59_2;
    int dim_rel_foi;
    int dim_S;
    int dim_S_1;
    int dim_S_2;
    int dim_S_preg;
    int dim_S_preg_1;
    int dim_S_preg_2;
    int dim_sample_rates;
    int dim_sample_transition_rates;
    int dim_T;
    int dim_T_1;
    int dim_T_2;
    int dim_T_preg;
    int dim_T_preg_1;
    int dim_T_preg_2;
    int dim_U;
    int dim_U_1;
    int dim_U_2;
    int dim_U_preg;
    int dim_U_preg_1;
    int dim_U_preg_2;
    int dim_x_I;
    int dim_Y;
    int dim_Y_1;
    int dim_Y_2;
    int dim_Y_preg;
    int dim_Y_preg_1;
    int dim_Y_preg_2;
    real_type DY;
    real_type eta;
    std::vector<real_type> fd;
    real_type fD0;
    std::vector<real_type> foi_age;
    std::vector<real_type> FOI_eq;
    real_type FOIv_eq;
    real_type ft;
    real_type fv;
    real_type gamma1;
    real_type gammaD;
    std::vector<real_type> het_wt;
    real_type IB0;
    real_type IC0;
    real_type ID0;
    std::vector<real_type> init_A;
    std::vector<real_type> init_A_preg;
    std::vector<real_type> init_D;
    std::vector<real_type> init_D_preg;
    real_type init_Ev;
    std::vector<real_type> init_FOI;
    std::vector<real_type> init_IB;
    std::vector<real_type> init_ICA;
    std::vector<real_type> init_ICM;
    std::vector<real_type> init_ID;
    real_type init_Iv;
    std::vector<real_type> init_P;
    std::vector<real_type> init_P_preg;
    std::vector<real_type> init_pregs;
    std::vector<real_type> init_S;
    std::vector<real_type> init_S_preg;
    real_type init_Sv;
    std::vector<real_type> init_T;
    std::vector<real_type> init_T_preg;
    std::vector<real_type> init_U;
    std::vector<real_type> init_U_preg;
    std::vector<real_type> initial_A;
    std::vector<real_type> initial_A_preg;
    real_type initial_betaa_td;
    std::vector<real_type> initial_D;
    std::vector<real_type> initial_D_preg;
    real_type initial_Ev;
    std::vector<real_type> initial_FOI;
    std::vector<real_type> initial_FOIv;
    std::vector<real_type> initial_IB;
    std::vector<real_type> initial_ICA;
    std::vector<real_type> initial_ICM;
    std::vector<real_type> initial_ID;
    std::vector<real_type> initial_ince_delay;
    real_type initial_Iv;
    std::vector<real_type> initial_P;
    std::vector<real_type> initial_P_preg;
    std::vector<real_type> initial_pregs;
    std::vector<real_type> initial_S;
    std::vector<real_type> initial_S_preg;
    real_type initial_Sv;
    std::vector<real_type> initial_T;
    std::vector<real_type> initial_T_preg;
    std::vector<real_type> initial_U;
    std::vector<real_type> initial_U_preg;
    real_type kB;
    real_type kC;
    real_type kD;
    int lag_rates;
    int lag_ratesMos;
    real_type mu;
    real_type mu0;
    real_type mv0;
    int na;
    int nh;
    int nrates;
    int offset_variable_A;
    int offset_variable_A_preg;
    int offset_variable_D;
    int offset_variable_D_preg;
    int offset_variable_FOI;
    int offset_variable_FOIv;
    int offset_variable_IB;
    int offset_variable_ICA;
    int offset_variable_ICM;
    int offset_variable_ID;
    int offset_variable_ince_delay;
    int offset_variable_P;
    int offset_variable_P_preg;
    int offset_variable_S;
    int offset_variable_S_preg;
    int offset_variable_T;
    int offset_variable_T_preg;
    int offset_variable_U;
    int offset_variable_U_preg;
    real_type omega;
    real_type p10;
    real_type p2;
    real_type phi0;
    real_type phi1;
    real_type PM;
    int preg_range;
    real_type rA;
    real_type rA_preg;
    real_type rD;
    std::vector<real_type> rel_foi;
    real_type rP;
    real_type rT;
    real_type rU;
    real_type rU_preg;
    std::vector<real_type> sample_rates;
    std::vector<real_type> sample_transition_rates;
    real_type surv;
    real_type tau1;
    real_type tau2;
    real_type uB;
    real_type uCA;
    real_type uD;
    real_type wane_rates;
    std::vector<real_type> x_I;
  };
  struct internal_type {
    std::vector<real_type> b;
    std::vector<real_type> cA;
    std::vector<real_type> clin_inc;
    std::vector<real_type> clin_inc_preg;
    std::vector<real_type> clin_inc0to5;
    std::vector<real_type> EIR;
    std::vector<real_type> FOI_lag;
    std::vector<real_type> FOIvijk;
    std::vector<real_type> IC;
    std::vector<real_type> init_ICM_pre;
    std::vector<real_type> p_det;
    std::vector<real_type> phi;
    std::vector<real_type> prev_cba;
    std::vector<real_type> prev_rdt;
    std::vector<real_type> prev0to59;
    std::vector<real_type> Y;
    std::vector<real_type> Y_preg;
  };
  mipodinmodel(const mode::pars_type<mipodinmodel>& pars) :
    shared(pars.shared), internal(pars.internal) {
  }
  size_t n_variables() {
    return shared->dim_A + shared->dim_A_preg + shared->dim_D + shared->dim_D_preg + shared->dim_FOI + shared->dim_FOIv + shared->dim_IB + shared->dim_ICA + shared->dim_ICM + shared->dim_ID + shared->dim_ince_delay + shared->dim_P + shared->dim_P_preg + shared->dim_pregs + shared->dim_S + shared->dim_S_preg + shared->dim_T + shared->dim_T_preg + shared->dim_U + shared->dim_U_preg + 4;
  }
  size_t n_output() {
    return 18;
  }
  std::vector<real_type> initial(size_t t) {
    std::vector<real_type> state(shared->dim_A + shared->dim_A_preg + shared->dim_D + shared->dim_D_preg + shared->dim_FOI + shared->dim_FOIv + shared->dim_IB + shared->dim_ICA + shared->dim_ICM + shared->dim_ID + shared->dim_ince_delay + shared->dim_P + shared->dim_P_preg + shared->dim_pregs + shared->dim_S + shared->dim_S_preg + shared->dim_T + shared->dim_T_preg + shared->dim_U + shared->dim_U_preg + 4);
    state[0] = shared->initial_betaa_td;
    state[1] = shared->initial_Sv;
    state[2] = shared->initial_Ev;
    state[3] = shared->initial_Iv;
    std::copy(shared->initial_pregs.begin(), shared->initial_pregs.end(), state.begin() + 4);
    std::copy(shared->initial_FOIv.begin(), shared->initial_FOIv.end(), state.begin() + shared->offset_variable_FOIv);
    std::copy(shared->initial_ince_delay.begin(), shared->initial_ince_delay.end(), state.begin() + shared->offset_variable_ince_delay);
    std::copy(shared->initial_S.begin(), shared->initial_S.end(), state.begin() + shared->offset_variable_S);
    std::copy(shared->initial_T.begin(), shared->initial_T.end(), state.begin() + shared->offset_variable_T);
    std::copy(shared->initial_D.begin(), shared->initial_D.end(), state.begin() + shared->offset_variable_D);
    std::copy(shared->initial_A.begin(), shared->initial_A.end(), state.begin() + shared->offset_variable_A);
    std::copy(shared->initial_U.begin(), shared->initial_U.end(), state.begin() + shared->offset_variable_U);
    std::copy(shared->initial_P.begin(), shared->initial_P.end(), state.begin() + shared->offset_variable_P);
    std::copy(shared->initial_S_preg.begin(), shared->initial_S_preg.end(), state.begin() + shared->offset_variable_S_preg);
    std::copy(shared->initial_T_preg.begin(), shared->initial_T_preg.end(), state.begin() + shared->offset_variable_T_preg);
    std::copy(shared->initial_D_preg.begin(), shared->initial_D_preg.end(), state.begin() + shared->offset_variable_D_preg);
    std::copy(shared->initial_A_preg.begin(), shared->initial_A_preg.end(), state.begin() + shared->offset_variable_A_preg);
    std::copy(shared->initial_U_preg.begin(), shared->initial_U_preg.end(), state.begin() + shared->offset_variable_U_preg);
    std::copy(shared->initial_P_preg.begin(), shared->initial_P_preg.end(), state.begin() + shared->offset_variable_P_preg);
    std::copy(shared->initial_ICM.begin(), shared->initial_ICM.end(), state.begin() + shared->offset_variable_ICM);
    std::copy(shared->initial_ICA.begin(), shared->initial_ICA.end(), state.begin() + shared->offset_variable_ICA);
    std::copy(shared->initial_IB.begin(), shared->initial_IB.end(), state.begin() + shared->offset_variable_IB);
    std::copy(shared->initial_ID.begin(), shared->initial_ID.end(), state.begin() + shared->offset_variable_ID);
    std::copy(shared->initial_FOI.begin(), shared->initial_FOI.end(), state.begin() + shared->offset_variable_FOI);
    return state;
  }
  void update_stochastic(double t, const std::vector<double>& state, rng_state_type& rng_state, std::vector<double>& state_next) {
    
  }
  void rhs(double t, const std::vector<double>& state, std::vector<double>& dstatedt) {
    const real_type * S = state.data() + shared->offset_variable_S;
    const real_type * T = state.data() + shared->offset_variable_T;
    const real_type * D = state.data() + shared->offset_variable_D;
    const real_type * A = state.data() + shared->offset_variable_A;
    const real_type * U = state.data() + shared->offset_variable_U;
    const real_type * P = state.data() + shared->offset_variable_P;
    const real_type * pregs = state.data() + 4;
    const real_type * S_preg = state.data() + shared->offset_variable_S_preg;
    const real_type * T_preg = state.data() + shared->offset_variable_T_preg;
    const real_type * D_preg = state.data() + shared->offset_variable_D_preg;
    const real_type * A_preg = state.data() + shared->offset_variable_A_preg;
    const real_type * U_preg = state.data() + shared->offset_variable_U_preg;
    const real_type * P_preg = state.data() + shared->offset_variable_P_preg;
    const real_type * ICM = state.data() + shared->offset_variable_ICM;
    const real_type * ICA = state.data() + shared->offset_variable_ICA;
    const real_type * IB = state.data() + shared->offset_variable_IB;
    const real_type * ID = state.data() + shared->offset_variable_ID;
    const real_type * FOI = state.data() + shared->offset_variable_FOI;
    const real_type * FOIv = state.data() + shared->offset_variable_FOIv;
    const real_type * ince_delay = state.data() + shared->offset_variable_ince_delay;
    const real_type betaa_td = state[0];
    const real_type Sv = state[1];
    const real_type Ev = state[2];
    const real_type Iv = state[3];
    dstatedt[0] = 0;
    real_type ince = FOIv[shared->lag_ratesMos - 1] * shared->lag_ratesMos / (real_type) shared->delayGam * Sv;
    {
       int i = 1;
       dstatedt[shared->offset_variable_ince_delay + i - 1] = ince - (shared->lag_ratesMos / (real_type) shared->delayMos) * ince_delay[0];
    }
    for (int i = 2; i <= shared->lag_ratesMos; ++i) {
      dstatedt[shared->offset_variable_ince_delay + i - 1] = (shared->lag_ratesMos / (real_type) shared->delayMos) * ince_delay[i - 1 - 1] - (shared->lag_ratesMos / (real_type) shared->delayMos) * ince_delay[i - 1];
    }
    real_type Ah = odin_sum2<real_type>(A, 0, shared->dim_A_1, 0, shared->dim_A_2, shared->dim_A_1);
    {
       int i = 1;
       for (int j = 1; j <= shared->nh; ++j) {
         dstatedt[shared->offset_variable_ICA + i - 1 + shared->dim_ICA_1 * (j - 1)] = (shared->lag_rates / (real_type) shared->dE) * FOI[shared->dim_FOI_12 * (shared->lag_rates - 1) + shared->dim_FOI_1 * (j - 1) + i - 1] / (real_type) ((shared->lag_rates / (real_type) shared->dE) * FOI[shared->dim_FOI_12 * (shared->lag_rates - 1) + shared->dim_FOI_1 * (j - 1) + i - 1] * shared->uCA + 1) - 1 / (real_type) shared->dCA * ICA[shared->dim_ICA_1 * (j - 1) + i - 1] - ICA[shared->dim_ICA_1 * (j - 1) + i - 1] / (real_type) shared->x_I[i - 1];
       }
    }
    for (int i = 2; i <= shared->na; ++i) {
      for (int j = 1; j <= shared->nh; ++j) {
        dstatedt[shared->offset_variable_ICA + i - 1 + shared->dim_ICA_1 * (j - 1)] = (shared->lag_rates / (real_type) shared->dE) * FOI[shared->dim_FOI_12 * (shared->lag_rates - 1) + shared->dim_FOI_1 * (j - 1) + i - 1] / (real_type) ((shared->lag_rates / (real_type) shared->dE) * FOI[shared->dim_FOI_12 * (shared->lag_rates - 1) + shared->dim_FOI_1 * (j - 1) + i - 1] * shared->uCA + 1) - 1 / (real_type) shared->dCA * ICA[shared->dim_ICA_1 * (j - 1) + i - 1] - (ICA[shared->dim_ICA_1 * (j - 1) + i - 1] - ICA[shared->dim_ICA_1 * (j - 1) + i - 1 - 1]) / (real_type) shared->x_I[i - 1];
      }
    }
    {
       int i = 1;
       for (int j = 1; j <= shared->nh; ++j) {
         dstatedt[shared->offset_variable_ID + i - 1 + shared->dim_ID_1 * (j - 1)] = (shared->lag_rates / (real_type) shared->dE) * FOI[shared->dim_FOI_12 * (shared->lag_rates - 1) + shared->dim_FOI_1 * (j - 1) + i - 1] / (real_type) ((shared->lag_rates / (real_type) shared->dE) * FOI[shared->dim_FOI_12 * (shared->lag_rates - 1) + shared->dim_FOI_1 * (j - 1) + i - 1] * shared->uD + 1) - ID[shared->dim_ID_1 * (j - 1) + i - 1] / (real_type) shared->dID - ID[shared->dim_ID_1 * (j - 1) + i - 1] / (real_type) shared->x_I[i - 1];
       }
    }
    for (int i = 2; i <= shared->na; ++i) {
      for (int j = 1; j <= shared->nh; ++j) {
        dstatedt[shared->offset_variable_ID + i - 1 + shared->dim_ID_1 * (j - 1)] = (shared->lag_rates / (real_type) shared->dE) * FOI[shared->dim_FOI_12 * (shared->lag_rates - 1) + shared->dim_FOI_1 * (j - 1) + i - 1] / (real_type) ((shared->lag_rates / (real_type) shared->dE) * FOI[shared->dim_FOI_12 * (shared->lag_rates - 1) + shared->dim_FOI_1 * (j - 1) + i - 1] * shared->uD + 1) - ID[shared->dim_ID_1 * (j - 1) + i - 1] / (real_type) shared->dID - (ID[shared->dim_ID_1 * (j - 1) + i - 1] - ID[shared->dim_ID_1 * (j - 1) + i - 1 - 1]) / (real_type) shared->x_I[i - 1];
      }
    }
    {
       int i = 1;
       for (int j = 1; j <= shared->nh; ++j) {
         dstatedt[shared->offset_variable_P + i - 1 + shared->dim_P_1 * (j - 1)] = shared->rT * T[shared->dim_T_1 * (j - 1) + i - 1] - shared->rP * P[shared->dim_P_1 * (j - 1) + i - 1] - (shared->eta + shared->age_rate[i - 1]) * P[shared->dim_P_1 * (j - 1) + i - 1];
       }
    }
    for (int i = 2; i <= shared->na; ++i) {
      for (int j = 1; j <= shared->nh; ++j) {
        dstatedt[shared->offset_variable_P + i - 1 + shared->dim_P_1 * (j - 1)] = shared->rT * T[shared->dim_T_1 * (j - 1) + i - 1] - shared->rP * P[shared->dim_P_1 * (j - 1) + i - 1] - (shared->eta + shared->age_rate[i - 1]) * P[shared->dim_P_1 * (j - 1) + i - 1] + shared->age_rate[i - 1 - 1] * P[shared->dim_P_1 * (j - 1) + i - 1 - 1];
      }
    }
    dstatedt[1] = - ince - shared->mu * Sv + betaa_td;
    {
       int i = 1;
       for (int j = 1; j <= shared->nh; ++j) {
         dstatedt[shared->offset_variable_U + i - 1 + shared->dim_U_1 * (j - 1)] = shared->rA * A[shared->dim_A_1 * (j - 1) + i - 1] - (shared->lag_rates / (real_type) shared->dE) * FOI[shared->dim_FOI_12 * (shared->lag_rates - 1) + shared->dim_FOI_1 * (j - 1) + i - 1] * U[shared->dim_U_1 * (j - 1) + i - 1] - shared->rU * U[shared->dim_U_1 * (j - 1) + i - 1] - (shared->eta + shared->age_rate[i - 1]) * U[shared->dim_U_1 * (j - 1) + i - 1];
       }
    }
    for (int i = 2; i <= shared->na; ++i) {
      for (int j = 1; j <= shared->nh; ++j) {
        dstatedt[shared->offset_variable_U + i - 1 + shared->dim_U_1 * (j - 1)] = shared->rA * A[shared->dim_A_1 * (j - 1) + i - 1] - (shared->lag_rates / (real_type) shared->dE) * FOI[shared->dim_FOI_12 * (shared->lag_rates - 1) + shared->dim_FOI_1 * (j - 1) + i - 1] * U[shared->dim_U_1 * (j - 1) + i - 1] - shared->rU * U[shared->dim_U_1 * (j - 1) + i - 1] - (shared->eta + shared->age_rate[i - 1]) * U[shared->dim_U_1 * (j - 1) + i - 1] + shared->age_rate[i - 1 - 1] * U[shared->dim_U_1 * (j - 1) + i - 1 - 1];
      }
    }
    real_type Dh = odin_sum2<real_type>(D, 0, shared->dim_D_1, 0, shared->dim_D_2, shared->dim_D_1);
    for (int i = 1; i <= shared->nh; ++i) {
      internal.init_ICM_pre[i - 1] = shared->PM * (ICA[shared->dim_ICA_1 * (i - 1) + shared->age20l - 1] + shared->age_20_factor * (ICA[shared->dim_ICA_1 * (i - 1) + shared->age20u - 1] - ICA[shared->dim_ICA_1 * (i - 1) + shared->age20l - 1]));
    }
    real_type Ph = odin_sum2<real_type>(P, 0, shared->dim_P_1, 0, shared->dim_P_2, shared->dim_P_1);
    real_type Sh = odin_sum2<real_type>(S, 0, shared->dim_S_1, 0, shared->dim_S_2, shared->dim_S_1);
    real_type Th = odin_sum2<real_type>(T, 0, shared->dim_T_1, 0, shared->dim_T_2, shared->dim_T_1);
    real_type Uh = odin_sum2<real_type>(U, 0, shared->dim_U_1, 0, shared->dim_U_2, shared->dim_U_1);
    for (int i = 1; i <= shared->na; ++i) {
      for (int j = 1; j <= shared->nh; ++j) {
        internal.b[i - 1 + shared->dim_b_1 * (j - 1)] = shared->b0 * ((1 - shared->b1) / (real_type) (1 + std::pow((IB[shared->dim_IB_1 * (j - 1) + i - 1] / (real_type) shared->IB0), shared->kB)) + shared->b1);
      }
    }
    {
       int i = 1;
       for (int j = 1; j <= shared->nh; ++j) {
         dstatedt[shared->offset_variable_ICM + i - 1 + shared->dim_ICM_1 * (j - 1)] = - 1 / (real_type) shared->dCM * ICM[shared->dim_ICM_1 * (j - 1) + i - 1] + (internal.init_ICM_pre[j - 1] - ICM[shared->dim_ICM_1 * (j - 1) + i - 1]) / (real_type) shared->x_I[i - 1];
       }
    }
    for (int i = 2; i <= shared->na; ++i) {
      for (int j = 1; j <= shared->nh; ++j) {
        dstatedt[shared->offset_variable_ICM + i - 1 + shared->dim_ICM_1 * (j - 1)] = - 1 / (real_type) shared->dCM * ICM[shared->dim_ICM_1 * (j - 1) + i - 1] - (ICM[shared->dim_ICM_1 * (j - 1) + i - 1] - ICM[shared->dim_ICM_1 * (j - 1) + i - 1 - 1]) / (real_type) shared->x_I[i - 1];
      }
    }
    for (int i = 1; i <= shared->preg_range; ++i) {
      for (int j = 1; j <= shared->nh; ++j) {
        dstatedt[shared->offset_variable_P_preg + i - 1 + shared->dim_P_preg_1 * (j - 1)] = shared->rT * T_preg[shared->dim_T_preg_1 * (j - 1) + i - 1] - shared->rP * P_preg[shared->dim_P_preg_1 * (j - 1) + i - 1];
      }
    }
    for (int i = 1; i <= shared->preg_range; ++i) {
      for (int j = 1; j <= shared->nh; ++j) {
        dstatedt[shared->offset_variable_S_preg + i - 1 + shared->dim_S_preg_1 * (j - 1)] = - (shared->lag_rates / (real_type) shared->dE) * FOI[shared->dim_FOI_12 * (shared->lag_rates - 1) + shared->dim_FOI_1 * (j - 1) + static_cast<int>((i + shared->agestart) - 1) - 1] * S_preg[shared->dim_S_preg_1 * (j - 1) + i - 1] + shared->rP * P_preg[shared->dim_P_preg_1 * (j - 1) + i - 1] + shared->rU_preg * U_preg[shared->dim_U_preg_1 * (j - 1) + i - 1];
      }
    }
    for (int i = 1; i <= shared->preg_range; ++i) {
      for (int j = 1; j <= shared->nh; ++j) {
        dstatedt[shared->offset_variable_U_preg + i - 1 + shared->dim_U_preg_1 * (j - 1)] = shared->rA_preg * A_preg[shared->dim_A_preg_1 * (j - 1) + i - 1] - (shared->lag_rates / (real_type) shared->dE) * FOI[shared->dim_FOI_12 * (shared->lag_rates - 1) + shared->dim_FOI_1 * (j - 1) + static_cast<int>((i + shared->agestart) - 1) - 1] * U_preg[shared->dim_U_preg_1 * (j - 1) + i - 1] - shared->rU_preg * U_preg[shared->dim_U_preg_1 * (j - 1) + i - 1];
      }
    }
    for (int i = 1; i <= shared->dim_EIR_1; ++i) {
      for (int j = 1; j <= shared->dim_EIR_2; ++j) {
        internal.EIR[i - 1 + shared->dim_EIR_1 * (j - 1)] = shared->rel_foi[j - 1] * shared->foi_age[i - 1] * Iv * shared->av0 / (real_type) shared->omega;
      }
    }
    real_type H = Sh + Th + Dh + Ah + Uh + Ph;
    for (int i = 1; i <= shared->dim_IC_1; ++i) {
      for (int j = 1; j <= shared->dim_IC_2; ++j) {
        internal.IC[i - 1 + shared->dim_IC_1 * (j - 1)] = ICM[shared->dim_ICM_1 * (j - 1) + i - 1] + ICA[shared->dim_ICA_1 * (j - 1) + i - 1];
      }
    }
    real_type incv = ince_delay[shared->lag_ratesMos - 1] * shared->lag_ratesMos / (real_type) shared->delayMos * shared->surv;
    for (int i = 1; i <= shared->dim_p_det_1; ++i) {
      for (int j = 1; j <= shared->dim_p_det_2; ++j) {
        internal.p_det[i - 1 + shared->dim_p_det_1 * (j - 1)] = shared->d1 + (1 - shared->d1) / (real_type) (1 + shared->fd[i - 1] * std::pow((ID[shared->dim_ID_1 * (j - 1) + i - 1] / (real_type) shared->ID0), shared->kD));
      }
    }
    for (int i = 1; i <= shared->na; ++i) {
      for (int j = 1; j <= shared->nh; ++j) {
        internal.Y[i - 1 + shared->dim_Y_1 * (j - 1)] = S[shared->dim_S_1 * (j - 1) + i - 1] + A[shared->dim_A_1 * (j - 1) + i - 1] + U[shared->dim_U_1 * (j - 1) + i - 1];
      }
    }
    for (int i = 1; i <= shared->dim_cA_1; ++i) {
      for (int j = 1; j <= shared->dim_cA_2; ++j) {
        internal.cA[i - 1 + shared->dim_cA_1 * (j - 1)] = shared->cU + (shared->cD - shared->cU) * std::pow(internal.p_det[shared->dim_p_det_1 * (j - 1) + i - 1], shared->gamma1);
      }
    }
    dstatedt[2] = ince - incv - shared->mu * Ev;
    {
       int i = 1;
       for (int j = 1; j <= shared->nh; ++j) {
         dstatedt[shared->offset_variable_IB + i - 1 + shared->dim_IB_1 * (j - 1)] = internal.EIR[shared->dim_EIR_1 * (j - 1) + i - 1] / (real_type) (internal.EIR[shared->dim_EIR_1 * (j - 1) + i - 1] * shared->uB + 1) - IB[shared->dim_IB_1 * (j - 1) + i - 1] / (real_type) shared->dB - IB[shared->dim_IB_1 * (j - 1) + i - 1] / (real_type) shared->x_I[i - 1];
       }
    }
    for (int i = 2; i <= shared->na; ++i) {
      for (int j = 1; j <= shared->nh; ++j) {
        dstatedt[shared->offset_variable_IB + i - 1 + shared->dim_IB_1 * (j - 1)] = internal.EIR[shared->dim_EIR_1 * (j - 1) + i - 1] / (real_type) (internal.EIR[shared->dim_EIR_1 * (j - 1) + i - 1] * shared->uB + 1) - IB[shared->dim_IB_1 * (j - 1) + i - 1] / (real_type) shared->dB - (IB[shared->dim_IB_1 * (j - 1) + i - 1] - IB[shared->dim_IB_1 * (j - 1) + i - 1 - 1]) / (real_type) shared->x_I[i - 1];
      }
    }
    dstatedt[3] = incv - shared->mu * Iv;
    {
       int i = 1;
       for (int j = 1; j <= shared->nh; ++j) {
         dstatedt[shared->offset_variable_S + i - 1 + shared->dim_S_1 * (j - 1)] = - (shared->lag_rates / (real_type) shared->dE) * FOI[shared->dim_FOI_12 * (shared->lag_rates - 1) + shared->dim_FOI_1 * (j - 1) + i - 1] * S[shared->dim_S_1 * (j - 1) + i - 1] + shared->rP * P[shared->dim_P_1 * (j - 1) + i - 1] + shared->rU * U[shared->dim_U_1 * (j - 1) + i - 1] + shared->eta * H * shared->het_wt[j - 1] - (shared->eta + shared->age_rate[i - 1]) * S[shared->dim_S_1 * (j - 1) + i - 1];
       }
    }
    for (int i = 2; i <= shared->na; ++i) {
      for (int j = 1; j <= shared->nh; ++j) {
        dstatedt[shared->offset_variable_S + i - 1 + shared->dim_S_1 * (j - 1)] = - (shared->lag_rates / (real_type) shared->dE) * FOI[shared->dim_FOI_12 * (shared->lag_rates - 1) + shared->dim_FOI_1 * (j - 1) + i - 1] * S[shared->dim_S_1 * (j - 1) + i - 1] + shared->rP * P[shared->dim_P_1 * (j - 1) + i - 1] + shared->rU * U[shared->dim_U_1 * (j - 1) + i - 1] - (shared->eta + shared->age_rate[i - 1]) * S[shared->dim_S_1 * (j - 1) + i - 1] + shared->age_rate[i - 1 - 1] * S[shared->dim_S_1 * (j - 1) + i - 1 - 1];
      }
    }
    for (int i = 1; i <= shared->na; ++i) {
      for (int j = 1; j <= shared->nh; ++j) {
        internal.FOI_lag[i - 1 + shared->dim_FOI_lag_1 * (j - 1)] = internal.EIR[shared->dim_EIR_1 * (j - 1) + i - 1] * ((IB[shared->dim_IB_1 * (j - 1) + i - 1] == 0 ? shared->b0 : internal.b[shared->dim_b_1 * (j - 1) + i - 1]));
      }
    }
    for (int i = 1; i <= shared->na; ++i) {
      for (int j = 1; j <= shared->nh; ++j) {
        internal.phi[i - 1 + shared->dim_phi_1 * (j - 1)] = shared->phi0 * ((1 - shared->phi1) / (real_type) (1 + std::pow((internal.IC[shared->dim_IC_1 * (j - 1) + i - 1] / (real_type) shared->IC0), shared->kC)) + shared->phi1);
      }
    }
    for (int i = 1; i <= shared->preg_range; ++i) {
      for (int j = 1; j <= shared->dim_prev_cba_2; ++j) {
        internal.prev_cba[i - 1 + shared->dim_prev_cba_1 * (j - 1)] = D_preg[shared->dim_D_preg_1 * (j - 1) + i - 1] + A_preg[shared->dim_A_preg_1 * (j - 1) + i - 1] * std::pow(internal.p_det[shared->dim_p_det_1 * (j - 1) + static_cast<int>((i + shared->agestart) - 1) - 1], shared->alphaA) + U_preg[shared->dim_U_preg_1 * (j - 1) + i - 1] * std::pow(internal.p_det[shared->dim_p_det_1 * (j - 1) + static_cast<int>((i + shared->agestart) - 1) - 1], shared->alphaU);
      }
    }
    for (int i = 1; i <= shared->preg_range; ++i) {
      for (int j = 1; j <= shared->nh; ++j) {
        internal.Y_preg[i - 1 + shared->dim_Y_preg_1 * (j - 1)] = S_preg[shared->dim_S_preg_1 * (j - 1) + i - 1] + A_preg[shared->dim_A_preg_1 * (j - 1) + i - 1] + U_preg[shared->dim_U_preg_1 * (j - 1) + i - 1];
      }
    }
    for (int i = 1; i <= shared->na; ++i) {
      for (int j = 1; j <= shared->nh; ++j) {
        internal.clin_inc[i - 1 + shared->dim_clin_inc_1 * (j - 1)] = internal.phi[shared->dim_phi_1 * (j - 1) + i - 1] * (shared->lag_rates / (real_type) shared->dE) * FOI[shared->dim_FOI_12 * (shared->lag_rates - 1) + shared->dim_FOI_1 * (j - 1) + i - 1] * internal.Y[shared->dim_Y_1 * (j - 1) + i - 1];
      }
    }
    for (int i = 1; i <= shared->preg_range; ++i) {
      for (int j = 1; j <= shared->nh; ++j) {
        internal.clin_inc_preg[i - 1 + shared->dim_clin_inc_preg_1 * (j - 1)] = internal.phi[shared->dim_phi_1 * (j - 1) + static_cast<int>((i + shared->agestart) - 1) - 1] * (shared->lag_rates / (real_type) shared->dE) * FOI[shared->dim_FOI_12 * (shared->lag_rates - 1) + shared->dim_FOI_1 * (j - 1) + static_cast<int>((i + shared->agestart) - 1) - 1] * internal.Y_preg[shared->dim_Y_preg_1 * (j - 1) + i - 1];
      }
    }
    {
       int i = 1;
       for (int j = 1; j <= shared->nh; ++j) {
         dstatedt[shared->offset_variable_A + i - 1 + shared->dim_A_1 * (j - 1)] = (1 - internal.phi[shared->dim_phi_1 * (j - 1) + i - 1]) * (shared->lag_rates / (real_type) shared->dE) * FOI[shared->dim_FOI_12 * (shared->lag_rates - 1) + shared->dim_FOI_1 * (j - 1) + i - 1] * internal.Y[shared->dim_Y_1 * (j - 1) + i - 1] - (shared->lag_rates / (real_type) shared->dE) * FOI[shared->dim_FOI_12 * (shared->lag_rates - 1) + shared->dim_FOI_1 * (j - 1) + i - 1] * A[shared->dim_A_1 * (j - 1) + i - 1] + shared->rD * D[shared->dim_D_1 * (j - 1) + i - 1] - shared->rA * A[shared->dim_A_1 * (j - 1) + i - 1] - (shared->eta + shared->age_rate[i - 1]) * A[shared->dim_A_1 * (j - 1) + i - 1];
       }
    }
    for (int i = 2; i <= shared->na; ++i) {
      for (int j = 1; j <= shared->nh; ++j) {
        dstatedt[shared->offset_variable_A + i - 1 + shared->dim_A_1 * (j - 1)] = (1 - internal.phi[shared->dim_phi_1 * (j - 1) + i - 1]) * (shared->lag_rates / (real_type) shared->dE) * FOI[shared->dim_FOI_12 * (shared->lag_rates - 1) + shared->dim_FOI_1 * (j - 1) + i - 1] * internal.Y[shared->dim_Y_1 * (j - 1) + i - 1] - (shared->lag_rates / (real_type) shared->dE) * FOI[shared->dim_FOI_12 * (shared->lag_rates - 1) + shared->dim_FOI_1 * (j - 1) + i - 1] * A[shared->dim_A_1 * (j - 1) + i - 1] + shared->rD * D[shared->dim_D_1 * (j - 1) + i - 1] - shared->rA * A[shared->dim_A_1 * (j - 1) + i - 1] - (shared->eta + shared->age_rate[i - 1]) * A[shared->dim_A_1 * (j - 1) + i - 1] + shared->age_rate[i - 1 - 1] * A[shared->dim_A_1 * (j - 1) + i - 1 - 1];
      }
    }
    for (int i = 1; i <= shared->preg_range; ++i) {
      for (int j = 1; j <= shared->nh; ++j) {
        dstatedt[shared->offset_variable_A_preg + i - 1 + shared->dim_A_preg_1 * (j - 1)] = (1 - internal.phi[shared->dim_phi_1 * (j - 1) + static_cast<int>((i + shared->agestart) - 1) - 1]) * (shared->lag_rates / (real_type) shared->dE) * FOI[shared->dim_FOI_12 * (shared->lag_rates - 1) + shared->dim_FOI_1 * (j - 1) + static_cast<int>((i + shared->agestart) - 1) - 1] * internal.Y_preg[shared->dim_Y_preg_1 * (j - 1) + i - 1] - (shared->lag_rates / (real_type) shared->dE) * FOI[shared->dim_FOI_12 * (shared->lag_rates - 1) + shared->dim_FOI_1 * (j - 1) + static_cast<int>((i + shared->agestart) - 1) - 1] * A_preg[shared->dim_A_preg_1 * (j - 1) + i - 1] + shared->rD * D_preg[shared->dim_D_preg_1 * (j - 1) + i - 1] - shared->rA_preg * A_preg[shared->dim_A_preg_1 * (j - 1) + i - 1];
      }
    }
    for (int i = 1; i <= shared->dim_FOI_1; ++i) {
      for (int j = 1; j <= shared->dim_FOI_2; ++j) {
        int k = 1;
        dstatedt[shared->offset_variable_FOI + i - 1 + shared->dim_FOI_1 * (j - 1) + shared->dim_FOI_12 * (k - 1)] = internal.FOI_lag[shared->dim_FOI_lag_1 * (j - 1) + i - 1] - (shared->lag_rates / (real_type) shared->dE) * FOI[shared->dim_FOI_12 * 0 + shared->dim_FOI_1 * (j - 1) + i - 1];
      }
    }
    for (int i = 1; i <= shared->dim_FOI_1; ++i) {
      for (int j = 1; j <= shared->dim_FOI_2; ++j) {
        for (int k = 2; k <= shared->lag_rates; ++k) {
          dstatedt[shared->offset_variable_FOI + i - 1 + shared->dim_FOI_1 * (j - 1) + shared->dim_FOI_12 * (k - 1)] = (shared->lag_rates / (real_type) shared->dE) * FOI[shared->dim_FOI_12 * (k - 1 - 1) + shared->dim_FOI_1 * (j - 1) + i - 1] - (shared->lag_rates / (real_type) shared->dE) * FOI[shared->dim_FOI_12 * (k - 1) + shared->dim_FOI_1 * (j - 1) + i - 1];
        }
      }
    }
    for (int i = 1; i <= shared->na; ++i) {
      for (int j = 1; j <= shared->nh; ++j) {
        internal.FOIvijk[i - 1 + shared->dim_FOIvijk_1 * (j - 1)] = shared->av0 * (shared->cT * T[shared->dim_T_1 * (j - 1) + i - 1] + shared->cD * D[shared->dim_D_1 * (j - 1) + i - 1] + internal.cA[shared->dim_cA_1 * (j - 1) + i - 1] * A[shared->dim_A_1 * (j - 1) + i - 1] + shared->cU * U[shared->dim_U_1 * (j - 1) + i - 1]) * shared->rel_foi[j - 1] * shared->foi_age[i - 1] / (real_type) shared->omega;
      }
    }
    real_type prev_pcr = odin_sum2<real_type>(internal.prev_cba.data(), 0, shared->dim_prev_cba_1, 0, shared->dim_prev_cba_2, shared->dim_prev_cba_1) / (real_type) odin_sum1<real_type>(shared->den.data(), shared->agestart - 1, shared->ageend);
    {
       int i = 1;
       for (int j = 1; j <= shared->nh; ++j) {
         dstatedt[shared->offset_variable_D + i - 1 + shared->dim_D_1 * (j - 1)] = (1 - shared->ft) * internal.clin_inc[shared->dim_clin_inc_1 * (j - 1) + i - 1] - shared->rD * D[shared->dim_D_1 * (j - 1) + i - 1] - (shared->eta + shared->age_rate[i - 1]) * D[shared->dim_D_1 * (j - 1) + i - 1];
       }
    }
    for (int i = 2; i <= shared->na; ++i) {
      for (int j = 1; j <= shared->nh; ++j) {
        dstatedt[shared->offset_variable_D + i - 1 + shared->dim_D_1 * (j - 1)] = (1 - shared->ft) * internal.clin_inc[shared->dim_clin_inc_1 * (j - 1) + i - 1] - shared->rD * D[shared->dim_D_1 * (j - 1) + i - 1] - (shared->eta + shared->age_rate[i - 1]) * D[shared->dim_D_1 * (j - 1) + i - 1] + shared->age_rate[i - 1 - 1] * D[shared->dim_D_1 * (j - 1) + i - 1 - 1];
      }
    }
    for (int i = 1; i <= shared->preg_range; ++i) {
      for (int j = 1; j <= shared->nh; ++j) {
        dstatedt[shared->offset_variable_D_preg + i - 1 + shared->dim_D_preg_1 * (j - 1)] = (1 - shared->ft) * internal.clin_inc_preg[shared->dim_clin_inc_preg_1 * (j - 1) + i - 1] - shared->rD * D_preg[shared->dim_D_preg_1 * (j - 1) + i - 1];
      }
    }
    {
       int i = 1;
       dstatedt[4 + i - 1] = prev_pcr * shared->sample_rates[i - 1] - (shared->sample_transition_rates[i - 1] + shared->wane_rates) * pregs[i - 1];
    }
    for (int i = 2; i <= shared->nrates; ++i) {
      dstatedt[4 + i - 1] = prev_pcr * shared->sample_rates[i - 1] + shared->sample_transition_rates[i - 1 - 1] * pregs[i - 1 - 1] - (shared->sample_transition_rates[i - 1] + shared->wane_rates) * pregs[i - 1];
    }
    {
       int i = 1;
       for (int j = 1; j <= shared->nh; ++j) {
         dstatedt[shared->offset_variable_T + i - 1 + shared->dim_T_1 * (j - 1)] = shared->ft * internal.clin_inc[shared->dim_clin_inc_1 * (j - 1) + i - 1] - shared->rT * T[shared->dim_T_1 * (j - 1) + i - 1] - (shared->eta + shared->age_rate[i - 1]) * T[shared->dim_T_1 * (j - 1) + i - 1];
       }
    }
    for (int i = 2; i <= shared->na; ++i) {
      for (int j = 1; j <= shared->nh; ++j) {
        dstatedt[shared->offset_variable_T + i - 1 + shared->dim_T_1 * (j - 1)] = shared->ft * internal.clin_inc[shared->dim_clin_inc_1 * (j - 1) + i - 1] - shared->rT * T[shared->dim_T_1 * (j - 1) + i - 1] - (shared->eta + shared->age_rate[i - 1]) * T[shared->dim_T_1 * (j - 1) + i - 1] + shared->age_rate[i - 1 - 1] * T[shared->dim_T_1 * (j - 1) + i - 1 - 1];
      }
    }
    for (int i = 1; i <= shared->preg_range; ++i) {
      for (int j = 1; j <= shared->nh; ++j) {
        dstatedt[shared->offset_variable_T_preg + i - 1 + shared->dim_T_preg_1 * (j - 1)] = shared->ft * internal.clin_inc_preg[shared->dim_clin_inc_preg_1 * (j - 1) + i - 1] - shared->rT * T_preg[shared->dim_T_preg_1 * (j - 1) + i - 1];
      }
    }
    real_type lag_FOIv = odin_sum1<real_type>(internal.FOIvijk.data(), 0, shared->dim_FOIvijk);
    {
       int i = 1;
       dstatedt[shared->offset_variable_FOIv + i - 1] = lag_FOIv - (shared->lag_ratesMos / (real_type) shared->delayGam) * FOIv[0];
    }
    for (int i = 2; i <= shared->lag_ratesMos; ++i) {
      dstatedt[shared->offset_variable_FOIv + i - 1] = (shared->lag_ratesMos / (real_type) shared->delayGam) * FOIv[i - 1 - 1] - (shared->lag_ratesMos / (real_type) shared->delayGam) * FOIv[i - 1];
    }
  }
  void output(double t, const std::vector<double>& state, std::vector<double>& output) {
    const real_type * S = state.data() + shared->offset_variable_S;
    const real_type * T = state.data() + shared->offset_variable_T;
    const real_type * D = state.data() + shared->offset_variable_D;
    const real_type * A = state.data() + shared->offset_variable_A;
    const real_type * U = state.data() + shared->offset_variable_U;
    const real_type * P = state.data() + shared->offset_variable_P;
    const real_type * S_preg = state.data() + shared->offset_variable_S_preg;
    const real_type * T_preg = state.data() + shared->offset_variable_T_preg;
    const real_type * D_preg = state.data() + shared->offset_variable_D_preg;
    const real_type * A_preg = state.data() + shared->offset_variable_A_preg;
    const real_type * U_preg = state.data() + shared->offset_variable_U_preg;
    const real_type * P_preg = state.data() + shared->offset_variable_P_preg;
    const real_type * ICM = state.data() + shared->offset_variable_ICM;
    const real_type * ICA = state.data() + shared->offset_variable_ICA;
    const real_type * ID = state.data() + shared->offset_variable_ID;
    const real_type * FOI = state.data() + shared->offset_variable_FOI;
    for (int i = 1; i <= shared->nh; ++i) {
      internal.init_ICM_pre[i - 1] = shared->PM * (ICA[shared->dim_ICA_1 * (i - 1) + shared->age20l - 1] + shared->age_20_factor * (ICA[shared->dim_ICA_1 * (i - 1) + shared->age20u - 1] - ICA[shared->dim_ICA_1 * (i - 1) + shared->age20l - 1]));
    }
    output[4] = odin_sum2<real_type>(A, 0, shared->dim_A_1, 0, shared->dim_A_2, shared->dim_A_1);
    output[3] = odin_sum2<real_type>(D, 0, shared->dim_D_1, 0, shared->dim_D_2, shared->dim_D_1);
    output[6] = odin_sum2<real_type>(P, 0, shared->dim_P_1, 0, shared->dim_P_2, shared->dim_P_1);
    output[1] = odin_sum2<real_type>(S, 0, shared->dim_S_1, 0, shared->dim_S_2, shared->dim_S_1);
    output[2] = odin_sum2<real_type>(T, 0, shared->dim_T_1, 0, shared->dim_T_2, shared->dim_T_1);
    output[5] = odin_sum2<real_type>(U, 0, shared->dim_U_1, 0, shared->dim_U_2, shared->dim_U_1);
    for (int i = 1; i <= shared->dim_IC_1; ++i) {
      for (int j = 1; j <= shared->dim_IC_2; ++j) {
        internal.IC[i - 1 + shared->dim_IC_1 * (j - 1)] = ICM[shared->dim_ICM_1 * (j - 1) + i - 1] + ICA[shared->dim_ICA_1 * (j - 1) + i - 1];
      }
    }
    output[10] = odin_sum2<real_type>(A_preg, 0, shared->dim_A_preg_1, 0, shared->dim_A_preg_2, shared->dim_A_preg_1);
    output[9] = odin_sum2<real_type>(D_preg, 0, shared->dim_D_preg_1, 0, shared->dim_D_preg_2, shared->dim_D_preg_1);
    output[0] = internal.init_ICM_pre[0];
    output[12] = odin_sum2<real_type>(P_preg, 0, shared->dim_P_preg_1, 0, shared->dim_P_preg_2, shared->dim_P_preg_1);
    output[7] = odin_sum2<real_type>(S_preg, 0, shared->dim_S_preg_1, 0, shared->dim_S_preg_2, shared->dim_S_preg_1);
    output[8] = odin_sum2<real_type>(T_preg, 0, shared->dim_T_preg_1, 0, shared->dim_T_preg_2, shared->dim_T_preg_1);
    output[11] = odin_sum2<real_type>(U_preg, 0, shared->dim_U_preg_1, 0, shared->dim_U_preg_2, shared->dim_U_preg_1);
    for (int i = 1; i <= shared->dim_p_det_1; ++i) {
      for (int j = 1; j <= shared->dim_p_det_2; ++j) {
        internal.p_det[i - 1 + shared->dim_p_det_1 * (j - 1)] = shared->d1 + (1 - shared->d1) / (real_type) (1 + shared->fd[i - 1] * std::pow((ID[shared->dim_ID_1 * (j - 1) + i - 1] / (real_type) shared->ID0), shared->kD));
      }
    }
    for (int i = 1; i <= shared->na; ++i) {
      for (int j = 1; j <= shared->nh; ++j) {
        internal.Y[i - 1 + shared->dim_Y_1 * (j - 1)] = S[shared->dim_S_1 * (j - 1) + i - 1] + A[shared->dim_A_1 * (j - 1) + i - 1] + U[shared->dim_U_1 * (j - 1) + i - 1];
      }
    }
    for (int i = 1; i <= shared->na; ++i) {
      for (int j = 1; j <= shared->nh; ++j) {
        internal.phi[i - 1 + shared->dim_phi_1 * (j - 1)] = shared->phi0 * ((1 - shared->phi1) / (real_type) (1 + std::pow((internal.IC[shared->dim_IC_1 * (j - 1) + i - 1] / (real_type) shared->IC0), shared->kC)) + shared->phi1);
      }
    }
    for (int i = 1; i <= shared->preg_range; ++i) {
      for (int j = 1; j <= shared->dim_prev_rdt_2; ++j) {
        internal.prev_rdt[i - 1 + shared->dim_prev_rdt_1 * (j - 1)] = D_preg[shared->dim_D_preg_1 * (j - 1) + i - 1] + A_preg[shared->dim_A_preg_1 * (j - 1) + i - 1] * internal.p_det[shared->dim_p_det_1 * (j - 1) + static_cast<int>((i + shared->agestart) - 1) - 1];
      }
    }
    for (int i = 1; i <= shared->age59; ++i) {
      for (int j = 1; j <= shared->dim_prev0to59_2; ++j) {
        internal.prev0to59[i - 1 + shared->dim_prev0to59_1 * (j - 1)] = T[shared->dim_T_1 * (j - 1) + i - 1] + D[shared->dim_D_1 * (j - 1) + i - 1] + A[shared->dim_A_1 * (j - 1) + i - 1] * internal.p_det[shared->dim_p_det_1 * (j - 1) + i - 1];
      }
    }
    for (int i = 1; i <= shared->preg_range; ++i) {
      for (int j = 1; j <= shared->nh; ++j) {
        internal.Y_preg[i - 1 + shared->dim_Y_preg_1 * (j - 1)] = S_preg[shared->dim_S_preg_1 * (j - 1) + i - 1] + A_preg[shared->dim_A_preg_1 * (j - 1) + i - 1] + U_preg[shared->dim_U_preg_1 * (j - 1) + i - 1];
      }
    }
    for (int i = 1; i <= shared->na; ++i) {
      for (int j = 1; j <= shared->nh; ++j) {
        internal.clin_inc[i - 1 + shared->dim_clin_inc_1 * (j - 1)] = internal.phi[shared->dim_phi_1 * (j - 1) + i - 1] * (shared->lag_rates / (real_type) shared->dE) * FOI[shared->dim_FOI_12 * (shared->lag_rates - 1) + shared->dim_FOI_1 * (j - 1) + i - 1] * internal.Y[shared->dim_Y_1 * (j - 1) + i - 1];
      }
    }
    output[14] = odin_sum2<real_type>(internal.prev0to59.data(), 0, shared->dim_prev0to59_1, 0, shared->dim_prev0to59_2, shared->dim_prev0to59_1) / (real_type) odin_sum1<real_type>(shared->den.data(), 0, shared->age59);
    output[17] = odin_sum2<real_type>(internal.prev_rdt.data(), 0, shared->dim_prev_rdt_1, 0, shared->dim_prev_rdt_2, shared->dim_prev_rdt_1) / (real_type) odin_sum1<real_type>(shared->den.data(), shared->agestart - 1, shared->ageend);
    output[13] = odin_sum2<real_type>(internal.Y_preg.data(), 0, shared->dim_Y_preg_1, 0, shared->dim_Y_preg_2, shared->dim_Y_preg_1);
    for (int i = 1; i <= shared->age05; ++i) {
      for (int j = 1; j <= shared->dim_clin_inc0to5_2; ++j) {
        internal.clin_inc0to5[i - 1 + shared->dim_clin_inc0to5_1 * (j - 1)] = internal.clin_inc[shared->dim_clin_inc_1 * (j - 1) + i - 1];
      }
    }
    output[16] = odin_sum2<real_type>(internal.clin_inc.data(), 0, shared->dim_clin_inc_1, 0, shared->dim_clin_inc_2, shared->dim_clin_inc_1);
    output[15] = odin_sum1<real_type>(internal.clin_inc0to5.data(), 0, shared->dim_clin_inc0to5) / (real_type) odin_sum1<real_type>(shared->den.data(), 0, shared->age05);
  }
private:
  std::shared_ptr<const shared_type> shared;
  internal_type internal;
};
template <typename real_type, typename container>
__host__ __device__ real_type odin_sum2(const container x, int from_i, int to_i, int from_j, int to_j, int dim_x_1) {
  real_type tot = 0.0;
  for (int j = from_j; j < to_j; ++j) {
    int jj = j * dim_x_1;
    for (int i = from_i; i < to_i; ++i) {
      tot += x[i + jj];
    }
  }
  return tot;
}
template <typename real_type, typename container>
__host__ __device__ real_type odin_sum3(const container x, int from_i, int to_i, int from_j, int to_j, int from_k, int to_k, int dim_x_1, int dim_x_12) {
  real_type tot = 0.0;
  for (int k = from_k; k < to_k; ++k) {
    int kk = k * dim_x_12;
    for (int j = from_j; j < to_j; ++j) {
      int jj = j * dim_x_1 + kk;
      for (int i = from_i; i < to_i; ++i) {
        tot += x[i + jj];
      }
    }
  }
  return tot;
}
#include <array>
#include <cpp11/R.hpp>
#include <cpp11/sexp.hpp>
#include <cpp11/doubles.hpp>
#include <cpp11/integers.hpp>
#include <cpp11/list.hpp>
#include <cpp11/strings.hpp>
#include <memory>
#include <vector>

template <typename T>
inline bool is_na(T x);

template <>
inline bool is_na(int x) {
  return x == NA_INTEGER;
}

template <>
inline bool is_na(double x) {
  return ISNA(x);
}

inline size_t object_length(cpp11::sexp x) {
  return ::Rf_xlength(x);
}

template <typename T>
void user_check_value(T value, const char *name, T min, T max) {
  if (is_na(value)) {
    cpp11::stop("'%s' must not be NA", name);
  }
  if (!is_na(min) && value < min) {
    cpp11::stop("Expected '%s' to be at least %g", name, (double) min);
  }
  if (!is_na(max) && value > max) {
    cpp11::stop("Expected '%s' to be at most %g", name, (double) max);
  }
}

template <typename T>
void user_check_array_value(const std::vector<T>& value, const char *name,
                            T min, T max) {
  for (auto& x : value) {
    user_check_value(x, name, min, max);
  }
}

inline size_t user_get_array_rank(cpp11::sexp x) {
  if (!::Rf_isArray(x)) {
    return 1;
  } else {
    cpp11::integers dim = cpp11::as_cpp<cpp11::integers>(x.attr("dim"));
    return dim.size();
  }
}

template <size_t N>
void user_check_array_rank(cpp11::sexp x, const char *name) {
  size_t rank = user_get_array_rank(x);
  if (rank != N) {
    if (N == 1) {
      cpp11::stop("Expected a vector for '%s'", name);
    } else if (N == 2) {
      cpp11::stop("Expected a matrix for '%s'", name);
    } else {
      cpp11::stop("Expected an array of rank %d for '%s'", N, name);
    }
  }
}

template <size_t N>
void user_check_array_dim(cpp11::sexp x, const char *name,
                          const std::array<int, N>& dim_expected) {
  cpp11::integers dim = cpp11::as_cpp<cpp11::integers>(x.attr("dim"));
  for (size_t i = 0; i < N; ++i) {
    if (dim[(int)i] != dim_expected[i]) {
      Rf_error("Incorrect size of dimension %d of '%s' (expected %d)",
               i + 1, name, dim_expected[i]);
    }
  }
}

template <>
inline void user_check_array_dim<1>(cpp11::sexp x, const char *name,
                                    const std::array<int, 1>& dim_expected) {
  if ((int)object_length(x) != dim_expected[0]) {
    cpp11::stop("Expected length %d value for '%s'", dim_expected[0], name);
  }
}

template <size_t N>
void user_set_array_dim(cpp11::sexp x, const char *name,
                        std::array<int, N>& dim) {
  cpp11::integers dim_given = cpp11::as_cpp<cpp11::integers>(x.attr("dim"));
  std::copy(dim_given.begin(), dim_given.end(), dim.begin());
}

template <>
inline void user_set_array_dim<1>(cpp11::sexp x, const char *name,
                                  std::array<int, 1>& dim) {
  dim[0] = object_length(x);
}

template <typename T>
T user_get_scalar(cpp11::list user, const char *name,
                  const T previous, T min, T max) {
  T ret = previous;
  cpp11::sexp x = user[name];
  if (x != R_NilValue) {
    if (object_length(x) != 1) {
      cpp11::stop("Expected a scalar numeric for '%s'", name);
    }
    // TODO: when we're getting out an integer this is a bit too relaxed
    if (TYPEOF(x) == REALSXP) {
      ret = cpp11::as_cpp<T>(x);
    } else if (TYPEOF(x) == INTSXP) {
      ret = cpp11::as_cpp<T>(x);
    } else {
      cpp11::stop("Expected a numeric value for %s", name);
    }
  }

  if (is_na(ret)) {
    cpp11::stop("Expected a value for '%s'", name);
  }
  user_check_value<T>(ret, name, min, max);
  return ret;
}

template <>
inline float user_get_scalar<float>(cpp11::list user, const char *name,
                                    const float previous, float min, float max) {
  double value = user_get_scalar<double>(user, name, previous, min, max);
  return static_cast<float>(value);
}

template <typename T>
std::vector<T> user_get_array_value(cpp11::sexp x, const char * name,
                                    T min, T max) {
  std::vector<T> ret = cpp11::as_cpp<std::vector<T>>(x);
  user_check_array_value<T>(ret, name, min, max);
  return ret;
}

template <typename T, size_t N>
std::vector<T> user_get_array_fixed(cpp11::list user, const char *name,
                                    const std::vector<T> previous,
                                    const std::array<int, N>& dim,
                                    T min, T max) {
  cpp11::sexp x = user[name];
  if (x == R_NilValue) {
    if (previous.size() == 0) {
      cpp11::stop("Expected a value for '%s'", name);
    }
    return previous;
  }

  user_check_array_rank<N>(x, name);
  user_check_array_dim<N>(x, name, dim);

  return user_get_array_value<T>(x, name, min, max);
}

template <typename T, size_t N>
std::vector<T> user_get_array_variable(cpp11::list user, const char *name,
                                       std::vector<T> previous,
                                       std::array<int, N>& dim,
                                       T min, T max) {
  cpp11::sexp x = user[name];
  if (x == R_NilValue) {
    if (previous.size() == 0) {
      cpp11::stop("Expected a value for '%s'", name);
    }
    return previous;
  }

  user_check_array_rank<N>(x, name);
  user_set_array_dim<N>(x, name, dim);

  return user_get_array_value<T>(x, name, min, max);
}

template <>
inline std::vector<float> user_get_array_value(cpp11::sexp x, const char * name,
                                               float min, float max) {
  // NOTE: possible under/overflow here for min/max because we've
  // downcast this.
  std::vector<double> value = user_get_array_value<double>(x, name, min, max);
  std::vector<float> ret(value.size());
  std::copy(value.begin(), value.end(), ret.begin());
  return ret;
}

// This is sum with inclusive "from", exclusive "to", following the
// same function in odin
template <typename real_type, typename container>
__host__ __device__
real_type odin_sum1(const container x, size_t from, size_t to) {
  real_type tot = 0.0;
  for (size_t i = from; i < to; ++i) {
    tot += x[i];
  }
  return tot;
}

inline cpp11::writable::integers integer_sequence(size_t from, size_t len) {
  cpp11::writable::integers ret(len);
  int* data = INTEGER(ret);
  for (size_t i = 0, j = from; i < len; ++i, ++j) {
    data[i] = j;
  }
  return ret;
}
namespace mode {
template<>
mode::pars_type<mipodinmodel> mode_pars<mipodinmodel>(cpp11::list user) {
  using real_type = typename mipodinmodel::real_type;
  auto shared = std::make_shared<mipodinmodel::shared_type>();
  mipodinmodel::internal_type internal;
  shared->initial_betaa_td = 0;
  shared->aD = NA_REAL;
  shared->age_20_factor = NA_REAL;
  shared->age05 = NA_INTEGER;
  shared->age20l = NA_INTEGER;
  shared->age20u = NA_INTEGER;
  shared->age59 = NA_INTEGER;
  shared->ageend = NA_REAL;
  shared->agestart = NA_REAL;
  shared->alphaA = NA_REAL;
  shared->alphaU = NA_REAL;
  shared->av0 = NA_REAL;
  shared->b0 = NA_REAL;
  shared->b1 = NA_REAL;
  shared->cD = NA_REAL;
  shared->cT = NA_REAL;
  shared->cU = NA_REAL;
  shared->d1 = NA_REAL;
  shared->dB = NA_REAL;
  shared->dCA = NA_REAL;
  shared->dCM = NA_REAL;
  shared->dE = NA_REAL;
  shared->delayGam = NA_REAL;
  shared->delayMos = NA_REAL;
  shared->dID = NA_REAL;
  shared->DY = NA_REAL;
  shared->eta = NA_REAL;
  shared->fD0 = NA_REAL;
  shared->FOIv_eq = NA_REAL;
  shared->ft = NA_REAL;
  shared->gamma1 = NA_REAL;
  shared->gammaD = NA_REAL;
  shared->IB0 = NA_REAL;
  shared->IC0 = NA_REAL;
  shared->ID0 = NA_REAL;
  shared->init_Ev = NA_REAL;
  shared->init_Iv = NA_REAL;
  shared->init_Sv = NA_REAL;
  shared->kB = NA_REAL;
  shared->kC = NA_REAL;
  shared->kD = NA_REAL;
  shared->lag_rates = NA_INTEGER;
  shared->lag_ratesMos = NA_INTEGER;
  shared->mu0 = NA_REAL;
  shared->mv0 = NA_REAL;
  shared->na = NA_INTEGER;
  shared->nh = NA_INTEGER;
  shared->nrates = NA_INTEGER;
  shared->omega = NA_REAL;
  shared->p10 = NA_REAL;
  shared->p2 = NA_REAL;
  shared->phi0 = NA_REAL;
  shared->phi1 = NA_REAL;
  shared->PM = NA_REAL;
  shared->rA = NA_REAL;
  shared->rA_preg = NA_REAL;
  shared->rD = NA_REAL;
  shared->rP = NA_REAL;
  shared->rT = NA_REAL;
  shared->rU = NA_REAL;
  shared->rU_preg = NA_REAL;
  shared->tau1 = NA_REAL;
  shared->tau2 = NA_REAL;
  shared->uB = NA_REAL;
  shared->uCA = NA_REAL;
  shared->uD = NA_REAL;
  shared->wane_rates = NA_REAL;
  shared->aD = user_get_scalar<real_type>(user, "aD", shared->aD, NA_REAL, NA_REAL);
  shared->age_20_factor = user_get_scalar<real_type>(user, "age_20_factor", shared->age_20_factor, NA_REAL, NA_REAL);
  shared->age05 = user_get_scalar<int>(user, "age05", shared->age05, NA_INTEGER, NA_INTEGER);
  shared->age20l = user_get_scalar<int>(user, "age20l", shared->age20l, NA_INTEGER, NA_INTEGER);
  shared->age20u = user_get_scalar<int>(user, "age20u", shared->age20u, NA_INTEGER, NA_INTEGER);
  shared->age59 = user_get_scalar<int>(user, "age59", shared->age59, NA_INTEGER, NA_INTEGER);
  shared->ageend = user_get_scalar<real_type>(user, "ageend", shared->ageend, NA_REAL, NA_REAL);
  shared->agestart = user_get_scalar<real_type>(user, "agestart", shared->agestart, NA_REAL, NA_REAL);
  shared->alphaA = user_get_scalar<real_type>(user, "alphaA", shared->alphaA, NA_REAL, NA_REAL);
  shared->alphaU = user_get_scalar<real_type>(user, "alphaU", shared->alphaU, NA_REAL, NA_REAL);
  shared->av0 = user_get_scalar<real_type>(user, "av0", shared->av0, NA_REAL, NA_REAL);
  shared->b0 = user_get_scalar<real_type>(user, "b0", shared->b0, NA_REAL, NA_REAL);
  shared->b1 = user_get_scalar<real_type>(user, "b1", shared->b1, NA_REAL, NA_REAL);
  shared->cD = user_get_scalar<real_type>(user, "cD", shared->cD, NA_REAL, NA_REAL);
  shared->cT = user_get_scalar<real_type>(user, "cT", shared->cT, NA_REAL, NA_REAL);
  shared->cU = user_get_scalar<real_type>(user, "cU", shared->cU, NA_REAL, NA_REAL);
  shared->d1 = user_get_scalar<real_type>(user, "d1", shared->d1, NA_REAL, NA_REAL);
  shared->dB = user_get_scalar<real_type>(user, "dB", shared->dB, NA_REAL, NA_REAL);
  shared->dCA = user_get_scalar<real_type>(user, "dCA", shared->dCA, NA_REAL, NA_REAL);
  shared->dCM = user_get_scalar<real_type>(user, "dCM", shared->dCM, NA_REAL, NA_REAL);
  shared->dE = user_get_scalar<real_type>(user, "dE", shared->dE, NA_REAL, NA_REAL);
  shared->delayGam = user_get_scalar<real_type>(user, "delayGam", shared->delayGam, NA_REAL, NA_REAL);
  shared->delayMos = user_get_scalar<real_type>(user, "delayMos", shared->delayMos, NA_REAL, NA_REAL);
  shared->dID = user_get_scalar<real_type>(user, "dID", shared->dID, NA_REAL, NA_REAL);
  shared->DY = user_get_scalar<real_type>(user, "DY", shared->DY, NA_REAL, NA_REAL);
  shared->eta = user_get_scalar<real_type>(user, "eta", shared->eta, NA_REAL, NA_REAL);
  shared->fD0 = user_get_scalar<real_type>(user, "fD0", shared->fD0, NA_REAL, NA_REAL);
  shared->FOIv_eq = user_get_scalar<real_type>(user, "FOIv_eq", shared->FOIv_eq, NA_REAL, NA_REAL);
  shared->ft = user_get_scalar<real_type>(user, "ft", shared->ft, NA_REAL, NA_REAL);
  shared->gamma1 = user_get_scalar<real_type>(user, "gamma1", shared->gamma1, NA_REAL, NA_REAL);
  shared->gammaD = user_get_scalar<real_type>(user, "gammaD", shared->gammaD, NA_REAL, NA_REAL);
  shared->IB0 = user_get_scalar<real_type>(user, "IB0", shared->IB0, NA_REAL, NA_REAL);
  shared->IC0 = user_get_scalar<real_type>(user, "IC0", shared->IC0, NA_REAL, NA_REAL);
  shared->ID0 = user_get_scalar<real_type>(user, "ID0", shared->ID0, NA_REAL, NA_REAL);
  shared->init_Ev = user_get_scalar<real_type>(user, "init_Ev", shared->init_Ev, NA_REAL, NA_REAL);
  shared->init_Iv = user_get_scalar<real_type>(user, "init_Iv", shared->init_Iv, NA_REAL, NA_REAL);
  shared->init_Sv = user_get_scalar<real_type>(user, "init_Sv", shared->init_Sv, NA_REAL, NA_REAL);
  shared->kB = user_get_scalar<real_type>(user, "kB", shared->kB, NA_REAL, NA_REAL);
  shared->kC = user_get_scalar<real_type>(user, "kC", shared->kC, NA_REAL, NA_REAL);
  shared->kD = user_get_scalar<real_type>(user, "kD", shared->kD, NA_REAL, NA_REAL);
  shared->lag_rates = user_get_scalar<int>(user, "lag_rates", shared->lag_rates, NA_INTEGER, NA_INTEGER);
  shared->lag_ratesMos = user_get_scalar<int>(user, "lag_ratesMos", shared->lag_ratesMos, NA_INTEGER, NA_INTEGER);
  shared->mu0 = user_get_scalar<real_type>(user, "mu0", shared->mu0, NA_REAL, NA_REAL);
  shared->mv0 = user_get_scalar<real_type>(user, "mv0", shared->mv0, NA_REAL, NA_REAL);
  shared->na = user_get_scalar<int>(user, "na", shared->na, NA_INTEGER, NA_INTEGER);
  shared->nh = user_get_scalar<int>(user, "nh", shared->nh, NA_INTEGER, NA_INTEGER);
  shared->nrates = user_get_scalar<int>(user, "nrates", shared->nrates, NA_INTEGER, NA_INTEGER);
  shared->omega = user_get_scalar<real_type>(user, "omega", shared->omega, NA_REAL, NA_REAL);
  shared->p10 = user_get_scalar<real_type>(user, "p10", shared->p10, NA_REAL, NA_REAL);
  shared->p2 = user_get_scalar<real_type>(user, "p2", shared->p2, NA_REAL, NA_REAL);
  shared->phi0 = user_get_scalar<real_type>(user, "phi0", shared->phi0, NA_REAL, NA_REAL);
  shared->phi1 = user_get_scalar<real_type>(user, "phi1", shared->phi1, NA_REAL, NA_REAL);
  shared->PM = user_get_scalar<real_type>(user, "PM", shared->PM, NA_REAL, NA_REAL);
  shared->rA = user_get_scalar<real_type>(user, "rA", shared->rA, NA_REAL, NA_REAL);
  shared->rA_preg = user_get_scalar<real_type>(user, "rA_preg", shared->rA_preg, NA_REAL, NA_REAL);
  shared->rD = user_get_scalar<real_type>(user, "rD", shared->rD, NA_REAL, NA_REAL);
  shared->rP = user_get_scalar<real_type>(user, "rP", shared->rP, NA_REAL, NA_REAL);
  shared->rT = user_get_scalar<real_type>(user, "rT", shared->rT, NA_REAL, NA_REAL);
  shared->rU = user_get_scalar<real_type>(user, "rU", shared->rU, NA_REAL, NA_REAL);
  shared->rU_preg = user_get_scalar<real_type>(user, "rU_preg", shared->rU_preg, NA_REAL, NA_REAL);
  shared->tau1 = user_get_scalar<real_type>(user, "tau1", shared->tau1, NA_REAL, NA_REAL);
  shared->tau2 = user_get_scalar<real_type>(user, "tau2", shared->tau2, NA_REAL, NA_REAL);
  shared->uB = user_get_scalar<real_type>(user, "uB", shared->uB, NA_REAL, NA_REAL);
  shared->uCA = user_get_scalar<real_type>(user, "uCA", shared->uCA, NA_REAL, NA_REAL);
  shared->uD = user_get_scalar<real_type>(user, "uD", shared->uD, NA_REAL, NA_REAL);
  shared->wane_rates = user_get_scalar<real_type>(user, "wane_rates", shared->wane_rates, NA_REAL, NA_REAL);
  shared->dim_A_1 = shared->na;
  shared->dim_A_2 = shared->nh;
  shared->dim_age = shared->na;
  shared->dim_age_rate = shared->na;
  shared->dim_b_1 = shared->na;
  shared->dim_b_2 = shared->nh;
  shared->dim_cA_1 = shared->na;
  shared->dim_cA_2 = shared->nh;
  shared->dim_clin_inc_1 = shared->na;
  shared->dim_clin_inc_2 = shared->nh;
  shared->dim_clin_inc0to5_1 = shared->age05;
  shared->dim_clin_inc0to5_2 = shared->nh;
  shared->dim_D_1 = shared->na;
  shared->dim_D_2 = shared->nh;
  shared->dim_den = shared->na;
  shared->dim_EIR_1 = shared->na;
  shared->dim_EIR_2 = shared->nh;
  shared->dim_fd = shared->na;
  shared->dim_FOI_1 = shared->na;
  shared->dim_FOI_2 = shared->nh;
  shared->dim_FOI_3 = shared->lag_rates;
  shared->dim_foi_age = shared->na;
  shared->dim_FOI_eq_1 = shared->na;
  shared->dim_FOI_eq_2 = shared->nh;
  shared->dim_FOI_lag_1 = shared->na;
  shared->dim_FOI_lag_2 = shared->nh;
  shared->dim_FOIv = shared->lag_ratesMos;
  shared->dim_FOIvijk_1 = shared->na;
  shared->dim_FOIvijk_2 = shared->nh;
  shared->dim_het_wt = shared->nh;
  shared->dim_IB_1 = shared->na;
  shared->dim_IB_2 = shared->nh;
  shared->dim_IC_1 = shared->na;
  shared->dim_IC_2 = shared->nh;
  shared->dim_ICA_1 = shared->na;
  shared->dim_ICA_2 = shared->nh;
  shared->dim_ICM_1 = shared->na;
  shared->dim_ICM_2 = shared->nh;
  shared->dim_ID_1 = shared->na;
  shared->dim_ID_2 = shared->nh;
  shared->dim_ince_delay = shared->lag_ratesMos;
  shared->dim_init_A_1 = shared->na;
  shared->dim_init_A_2 = shared->nh;
  shared->dim_init_D_1 = shared->na;
  shared->dim_init_D_2 = shared->nh;
  shared->dim_init_FOI_1 = shared->na;
  shared->dim_init_FOI_2 = shared->nh;
  shared->dim_init_FOI_3 = shared->lag_rates;
  shared->dim_init_IB_1 = shared->na;
  shared->dim_init_IB_2 = shared->nh;
  shared->dim_init_ICA_1 = shared->na;
  shared->dim_init_ICA_2 = shared->nh;
  shared->dim_init_ICM_1 = shared->na;
  shared->dim_init_ICM_2 = shared->nh;
  shared->dim_init_ICM_pre = shared->nh;
  shared->dim_init_ID_1 = shared->na;
  shared->dim_init_ID_2 = shared->nh;
  shared->dim_init_P_1 = shared->na;
  shared->dim_init_P_2 = shared->nh;
  shared->dim_init_pregs = shared->nrates;
  shared->dim_init_S_1 = shared->na;
  shared->dim_init_S_2 = shared->nh;
  shared->dim_init_T_1 = shared->na;
  shared->dim_init_T_2 = shared->nh;
  shared->dim_init_U_1 = shared->na;
  shared->dim_init_U_2 = shared->nh;
  shared->dim_P_1 = shared->na;
  shared->dim_P_2 = shared->nh;
  shared->dim_p_det_1 = shared->na;
  shared->dim_p_det_2 = shared->nh;
  shared->dim_phi_1 = shared->na;
  shared->dim_phi_2 = shared->nh;
  shared->dim_pregs = shared->nrates;
  shared->dim_prev0to59_1 = shared->age59;
  shared->dim_prev0to59_2 = shared->nh;
  shared->dim_rel_foi = shared->nh;
  shared->dim_S_1 = shared->na;
  shared->dim_S_2 = shared->nh;
  shared->dim_sample_rates = shared->nrates;
  shared->dim_sample_transition_rates = shared->nrates;
  shared->dim_T_1 = shared->na;
  shared->dim_T_2 = shared->nh;
  shared->dim_U_1 = shared->na;
  shared->dim_U_2 = shared->nh;
  shared->dim_x_I = shared->na;
  shared->dim_Y_1 = shared->na;
  shared->dim_Y_2 = shared->nh;
  shared->fv = 1 / (real_type) (shared->tau1 + shared->tau2);
  shared->initial_Ev = shared->init_Ev * shared->mv0;
  shared->initial_Iv = shared->init_Iv * shared->mv0;
  shared->initial_Sv = shared->init_Sv * shared->mv0;
  shared->preg_range = (shared->ageend - shared->agestart) + 1;
  shared->age = user_get_array_fixed<real_type, 1>(user, "age", shared->age, {shared->dim_age}, NA_REAL, NA_REAL);
  shared->age_rate = user_get_array_fixed<real_type, 1>(user, "age_rate", shared->age_rate, {shared->dim_age_rate}, NA_REAL, NA_REAL);
  shared->fd = std::vector<real_type>(shared->dim_fd);
  internal.init_ICM_pre = std::vector<real_type>(shared->dim_init_ICM_pre);
  shared->initial_FOIv = std::vector<real_type>(shared->dim_FOIv);
  shared->initial_ince_delay = std::vector<real_type>(shared->dim_ince_delay);
  shared->initial_pregs = std::vector<real_type>(shared->dim_pregs);
  shared->den = user_get_array_fixed<real_type, 1>(user, "den", shared->den, {shared->dim_den}, NA_REAL, NA_REAL);
  shared->dim_A = shared->dim_A_1 * shared->dim_A_2;
  shared->dim_A_preg_1 = shared->preg_range;
  shared->dim_A_preg_2 = shared->nh;
  shared->dim_b = shared->dim_b_1 * shared->dim_b_2;
  shared->dim_cA = shared->dim_cA_1 * shared->dim_cA_2;
  shared->dim_clin_inc = shared->dim_clin_inc_1 * shared->dim_clin_inc_2;
  shared->dim_clin_inc_preg_1 = shared->preg_range;
  shared->dim_clin_inc_preg_2 = shared->nh;
  shared->dim_clin_inc0to5 = shared->dim_clin_inc0to5_1 * shared->dim_clin_inc0to5_2;
  shared->dim_D = shared->dim_D_1 * shared->dim_D_2;
  shared->dim_D_preg_1 = shared->preg_range;
  shared->dim_D_preg_2 = shared->nh;
  shared->dim_EIR = shared->dim_EIR_1 * shared->dim_EIR_2;
  shared->dim_FOI = shared->dim_FOI_1 * shared->dim_FOI_2 * shared->dim_FOI_3;
  shared->dim_FOI_12 = shared->dim_FOI_1 * shared->dim_FOI_2;
  shared->dim_FOI_eq = shared->dim_FOI_eq_1 * shared->dim_FOI_eq_2;
  shared->dim_FOI_lag = shared->dim_FOI_lag_1 * shared->dim_FOI_lag_2;
  shared->dim_FOIvijk = shared->dim_FOIvijk_1 * shared->dim_FOIvijk_2;
  shared->dim_IB = shared->dim_IB_1 * shared->dim_IB_2;
  shared->dim_IC = shared->dim_IC_1 * shared->dim_IC_2;
  shared->dim_ICA = shared->dim_ICA_1 * shared->dim_ICA_2;
  shared->dim_ICM = shared->dim_ICM_1 * shared->dim_ICM_2;
  shared->dim_ID = shared->dim_ID_1 * shared->dim_ID_2;
  shared->dim_init_A = shared->dim_init_A_1 * shared->dim_init_A_2;
  shared->dim_init_A_preg_1 = shared->preg_range;
  shared->dim_init_A_preg_2 = shared->nh;
  shared->dim_init_D = shared->dim_init_D_1 * shared->dim_init_D_2;
  shared->dim_init_D_preg_1 = shared->preg_range;
  shared->dim_init_D_preg_2 = shared->nh;
  shared->dim_init_FOI = shared->dim_init_FOI_1 * shared->dim_init_FOI_2 * shared->dim_init_FOI_3;
  shared->dim_init_FOI_12 = shared->dim_init_FOI_1 * shared->dim_init_FOI_2;
  shared->dim_init_IB = shared->dim_init_IB_1 * shared->dim_init_IB_2;
  shared->dim_init_ICA = shared->dim_init_ICA_1 * shared->dim_init_ICA_2;
  shared->dim_init_ICM = shared->dim_init_ICM_1 * shared->dim_init_ICM_2;
  shared->dim_init_ID = shared->dim_init_ID_1 * shared->dim_init_ID_2;
  shared->dim_init_P = shared->dim_init_P_1 * shared->dim_init_P_2;
  shared->dim_init_P_preg_1 = shared->preg_range;
  shared->dim_init_P_preg_2 = shared->nh;
  shared->dim_init_S = shared->dim_init_S_1 * shared->dim_init_S_2;
  shared->dim_init_S_preg_1 = shared->preg_range;
  shared->dim_init_S_preg_2 = shared->nh;
  shared->dim_init_T = shared->dim_init_T_1 * shared->dim_init_T_2;
  shared->dim_init_T_preg_1 = shared->preg_range;
  shared->dim_init_T_preg_2 = shared->nh;
  shared->dim_init_U = shared->dim_init_U_1 * shared->dim_init_U_2;
  shared->dim_init_U_preg_1 = shared->preg_range;
  shared->dim_init_U_preg_2 = shared->nh;
  shared->dim_P = shared->dim_P_1 * shared->dim_P_2;
  shared->dim_p_det = shared->dim_p_det_1 * shared->dim_p_det_2;
  shared->dim_P_preg_1 = shared->preg_range;
  shared->dim_P_preg_2 = shared->nh;
  shared->dim_phi = shared->dim_phi_1 * shared->dim_phi_2;
  shared->dim_prev_cba_1 = shared->preg_range;
  shared->dim_prev_cba_2 = shared->nh;
  shared->dim_prev_rdt_1 = shared->preg_range;
  shared->dim_prev_rdt_2 = shared->nh;
  shared->dim_prev0to59 = shared->dim_prev0to59_1 * shared->dim_prev0to59_2;
  shared->dim_S = shared->dim_S_1 * shared->dim_S_2;
  shared->dim_S_preg_1 = shared->preg_range;
  shared->dim_S_preg_2 = shared->nh;
  shared->dim_T = shared->dim_T_1 * shared->dim_T_2;
  shared->dim_T_preg_1 = shared->preg_range;
  shared->dim_T_preg_2 = shared->nh;
  shared->dim_U = shared->dim_U_1 * shared->dim_U_2;
  shared->dim_U_preg_1 = shared->preg_range;
  shared->dim_U_preg_2 = shared->nh;
  shared->dim_Y = shared->dim_Y_1 * shared->dim_Y_2;
  shared->dim_Y_preg_1 = shared->preg_range;
  shared->dim_Y_preg_2 = shared->nh;
  shared->foi_age = user_get_array_fixed<real_type, 1>(user, "foi_age", shared->foi_age, {shared->dim_foi_age}, NA_REAL, NA_REAL);
  shared->het_wt = user_get_array_fixed<real_type, 1>(user, "het_wt", shared->het_wt, {shared->dim_het_wt}, NA_REAL, NA_REAL);
  shared->init_pregs = user_get_array_fixed<real_type, 1>(user, "init_pregs", shared->init_pregs, {shared->dim_init_pregs}, NA_REAL, NA_REAL);
  for (int i = 1; i <= shared->dim_FOIv; ++i) {
    shared->initial_FOIv[i - 1] = shared->FOIv_eq * shared->delayGam / (real_type) shared->lag_ratesMos;
  }
  for (int i = 1; i <= shared->dim_ince_delay; ++i) {
    shared->initial_ince_delay[i - 1] = shared->FOIv_eq * shared->init_Sv * shared->mv0 * shared->delayMos / (real_type) shared->lag_ratesMos;
  }
  shared->mu = - shared->fv * std::log(shared->p10 * shared->p2);
  shared->offset_variable_FOIv = shared->dim_pregs + 4;
  shared->offset_variable_ince_delay = shared->dim_FOIv + shared->dim_pregs + 4;
  shared->offset_variable_S = shared->dim_FOIv + shared->dim_ince_delay + shared->dim_pregs + 4;
  shared->rel_foi = user_get_array_fixed<real_type, 1>(user, "rel_foi", shared->rel_foi, {shared->dim_rel_foi}, NA_REAL, NA_REAL);
  shared->sample_rates = user_get_array_fixed<real_type, 1>(user, "sample_rates", shared->sample_rates, {shared->dim_sample_rates}, NA_REAL, NA_REAL);
  shared->sample_transition_rates = user_get_array_fixed<real_type, 1>(user, "sample_transition_rates", shared->sample_transition_rates, {shared->dim_sample_transition_rates}, NA_REAL, NA_REAL);
  shared->x_I = user_get_array_fixed<real_type, 1>(user, "x_I", shared->x_I, {shared->dim_x_I}, NA_REAL, NA_REAL);
  internal.b = std::vector<real_type>(shared->dim_b);
  internal.cA = std::vector<real_type>(shared->dim_cA);
  internal.clin_inc = std::vector<real_type>(shared->dim_clin_inc);
  internal.clin_inc0to5 = std::vector<real_type>(shared->dim_clin_inc0to5);
  internal.EIR = std::vector<real_type>(shared->dim_EIR);
  internal.FOI_lag = std::vector<real_type>(shared->dim_FOI_lag);
  internal.FOIvijk = std::vector<real_type>(shared->dim_FOIvijk);
  internal.IC = std::vector<real_type>(shared->dim_IC);
  shared->init_FOI = std::vector<real_type>(shared->dim_init_FOI);
  shared->initial_A = std::vector<real_type>(shared->dim_A);
  shared->initial_D = std::vector<real_type>(shared->dim_D);
  shared->initial_FOI = std::vector<real_type>(shared->dim_FOI);
  shared->initial_IB = std::vector<real_type>(shared->dim_IB);
  shared->initial_ICA = std::vector<real_type>(shared->dim_ICA);
  shared->initial_ICM = std::vector<real_type>(shared->dim_ICM);
  shared->initial_ID = std::vector<real_type>(shared->dim_ID);
  shared->initial_P = std::vector<real_type>(shared->dim_P);
  shared->initial_S = std::vector<real_type>(shared->dim_S);
  shared->initial_T = std::vector<real_type>(shared->dim_T);
  shared->initial_U = std::vector<real_type>(shared->dim_U);
  internal.p_det = std::vector<real_type>(shared->dim_p_det);
  internal.phi = std::vector<real_type>(shared->dim_phi);
  internal.prev0to59 = std::vector<real_type>(shared->dim_prev0to59);
  internal.Y = std::vector<real_type>(shared->dim_Y);
  shared->dim_A_preg = shared->dim_A_preg_1 * shared->dim_A_preg_2;
  shared->dim_clin_inc_preg = shared->dim_clin_inc_preg_1 * shared->dim_clin_inc_preg_2;
  shared->dim_D_preg = shared->dim_D_preg_1 * shared->dim_D_preg_2;
  shared->dim_init_A_preg = shared->dim_init_A_preg_1 * shared->dim_init_A_preg_2;
  shared->dim_init_D_preg = shared->dim_init_D_preg_1 * shared->dim_init_D_preg_2;
  shared->dim_init_P_preg = shared->dim_init_P_preg_1 * shared->dim_init_P_preg_2;
  shared->dim_init_S_preg = shared->dim_init_S_preg_1 * shared->dim_init_S_preg_2;
  shared->dim_init_T_preg = shared->dim_init_T_preg_1 * shared->dim_init_T_preg_2;
  shared->dim_init_U_preg = shared->dim_init_U_preg_1 * shared->dim_init_U_preg_2;
  shared->dim_P_preg = shared->dim_P_preg_1 * shared->dim_P_preg_2;
  shared->dim_prev_cba = shared->dim_prev_cba_1 * shared->dim_prev_cba_2;
  shared->dim_prev_rdt = shared->dim_prev_rdt_1 * shared->dim_prev_rdt_2;
  shared->dim_S_preg = shared->dim_S_preg_1 * shared->dim_S_preg_2;
  shared->dim_T_preg = shared->dim_T_preg_1 * shared->dim_T_preg_2;
  shared->dim_U_preg = shared->dim_U_preg_1 * shared->dim_U_preg_2;
  shared->dim_Y_preg = shared->dim_Y_preg_1 * shared->dim_Y_preg_2;
  for (int i = 1; i <= shared->na; ++i) {
    shared->fd[i - 1] = 1 - (1 - shared->fD0) / (real_type) (1 + std::pow((shared->age[i - 1] / (real_type) shared->aD), shared->gammaD));
  }
  shared->FOI_eq = user_get_array_fixed<real_type, 2>(user, "FOI_eq", shared->FOI_eq, {shared->dim_FOI_eq_1, shared->dim_FOI_eq_2}, NA_REAL, NA_REAL);
  shared->init_A = user_get_array_fixed<real_type, 2>(user, "init_A", shared->init_A, {shared->dim_init_A_1, shared->dim_init_A_2}, NA_REAL, NA_REAL);
  shared->init_D = user_get_array_fixed<real_type, 2>(user, "init_D", shared->init_D, {shared->dim_init_D_1, shared->dim_init_D_2}, NA_REAL, NA_REAL);
  shared->init_IB = user_get_array_fixed<real_type, 2>(user, "init_IB", shared->init_IB, {shared->dim_init_IB_1, shared->dim_init_IB_2}, NA_REAL, NA_REAL);
  shared->init_ICA = user_get_array_fixed<real_type, 2>(user, "init_ICA", shared->init_ICA, {shared->dim_init_ICA_1, shared->dim_init_ICA_2}, NA_REAL, NA_REAL);
  shared->init_ICM = user_get_array_fixed<real_type, 2>(user, "init_ICM", shared->init_ICM, {shared->dim_init_ICM_1, shared->dim_init_ICM_2}, NA_REAL, NA_REAL);
  shared->init_ID = user_get_array_fixed<real_type, 2>(user, "init_ID", shared->init_ID, {shared->dim_init_ID_1, shared->dim_init_ID_2}, NA_REAL, NA_REAL);
  shared->init_P = user_get_array_fixed<real_type, 2>(user, "init_P", shared->init_P, {shared->dim_init_P_1, shared->dim_init_P_2}, NA_REAL, NA_REAL);
  shared->init_S = user_get_array_fixed<real_type, 2>(user, "init_S", shared->init_S, {shared->dim_init_S_1, shared->dim_init_S_2}, NA_REAL, NA_REAL);
  shared->init_T = user_get_array_fixed<real_type, 2>(user, "init_T", shared->init_T, {shared->dim_init_T_1, shared->dim_init_T_2}, NA_REAL, NA_REAL);
  shared->init_U = user_get_array_fixed<real_type, 2>(user, "init_U", shared->init_U, {shared->dim_init_U_1, shared->dim_init_U_2}, NA_REAL, NA_REAL);
  for (int i = 1; i <= shared->dim_pregs; ++i) {
    shared->initial_pregs[i - 1] = shared->init_pregs[i - 1];
  }
  shared->offset_variable_A = shared->dim_D + shared->dim_FOIv + shared->dim_ince_delay + shared->dim_pregs + shared->dim_S + shared->dim_T + 4;
  shared->offset_variable_D = shared->dim_FOIv + shared->dim_ince_delay + shared->dim_pregs + shared->dim_S + shared->dim_T + 4;
  shared->offset_variable_P = shared->dim_A + shared->dim_D + shared->dim_FOIv + shared->dim_ince_delay + shared->dim_pregs + shared->dim_S + shared->dim_T + shared->dim_U + 4;
  shared->offset_variable_S_preg = shared->dim_A + shared->dim_D + shared->dim_FOIv + shared->dim_ince_delay + shared->dim_P + shared->dim_pregs + shared->dim_S + shared->dim_T + shared->dim_U + 4;
  shared->offset_variable_T = shared->dim_FOIv + shared->dim_ince_delay + shared->dim_pregs + shared->dim_S + 4;
  shared->offset_variable_U = shared->dim_A + shared->dim_D + shared->dim_FOIv + shared->dim_ince_delay + shared->dim_pregs + shared->dim_S + shared->dim_T + 4;
  shared->surv = std::exp(- shared->mu * shared->delayMos);
  internal.clin_inc_preg = std::vector<real_type>(shared->dim_clin_inc_preg);
  shared->initial_A_preg = std::vector<real_type>(shared->dim_A_preg);
  shared->initial_D_preg = std::vector<real_type>(shared->dim_D_preg);
  shared->initial_P_preg = std::vector<real_type>(shared->dim_P_preg);
  shared->initial_S_preg = std::vector<real_type>(shared->dim_S_preg);
  shared->initial_T_preg = std::vector<real_type>(shared->dim_T_preg);
  shared->initial_U_preg = std::vector<real_type>(shared->dim_U_preg);
  internal.prev_cba = std::vector<real_type>(shared->dim_prev_cba);
  internal.prev_rdt = std::vector<real_type>(shared->dim_prev_rdt);
  internal.Y_preg = std::vector<real_type>(shared->dim_Y_preg);
  shared->init_A_preg = user_get_array_fixed<real_type, 2>(user, "init_A_preg", shared->init_A_preg, {shared->dim_init_A_preg_1, shared->dim_init_A_preg_2}, NA_REAL, NA_REAL);
  shared->init_D_preg = user_get_array_fixed<real_type, 2>(user, "init_D_preg", shared->init_D_preg, {shared->dim_init_D_preg_1, shared->dim_init_D_preg_2}, NA_REAL, NA_REAL);
  for (int i = 1; i <= shared->dim_init_FOI_1; ++i) {
    for (int j = 1; j <= shared->dim_init_FOI_2; ++j) {
      for (int k = 1; k <= shared->dim_init_FOI_3; ++k) {
        shared->init_FOI[i - 1 + shared->dim_init_FOI_1 * (j - 1) + shared->dim_init_FOI_12 * (k - 1)] = shared->FOI_eq[shared->dim_FOI_eq_1 * (j - 1) + i - 1] * shared->dE / (real_type) shared->lag_rates;
      }
    }
  }
  shared->init_P_preg = user_get_array_fixed<real_type, 2>(user, "init_P_preg", shared->init_P_preg, {shared->dim_init_P_preg_1, shared->dim_init_P_preg_2}, NA_REAL, NA_REAL);
  shared->init_S_preg = user_get_array_fixed<real_type, 2>(user, "init_S_preg", shared->init_S_preg, {shared->dim_init_S_preg_1, shared->dim_init_S_preg_2}, NA_REAL, NA_REAL);
  shared->init_T_preg = user_get_array_fixed<real_type, 2>(user, "init_T_preg", shared->init_T_preg, {shared->dim_init_T_preg_1, shared->dim_init_T_preg_2}, NA_REAL, NA_REAL);
  shared->init_U_preg = user_get_array_fixed<real_type, 2>(user, "init_U_preg", shared->init_U_preg, {shared->dim_init_U_preg_1, shared->dim_init_U_preg_2}, NA_REAL, NA_REAL);
  for (int i = 1; i <= shared->dim_A_1; ++i) {
    for (int j = 1; j <= shared->dim_A_2; ++j) {
      shared->initial_A[i - 1 + shared->dim_A_1 * (j - 1)] = shared->init_A[shared->dim_init_A_1 * (j - 1) + i - 1];
    }
  }
  for (int i = 1; i <= shared->dim_D_1; ++i) {
    for (int j = 1; j <= shared->dim_D_2; ++j) {
      shared->initial_D[i - 1 + shared->dim_D_1 * (j - 1)] = shared->init_D[shared->dim_init_D_1 * (j - 1) + i - 1];
    }
  }
  for (int i = 1; i <= shared->dim_IB_1; ++i) {
    for (int j = 1; j <= shared->dim_IB_2; ++j) {
      shared->initial_IB[i - 1 + shared->dim_IB_1 * (j - 1)] = shared->init_IB[shared->dim_init_IB_1 * (j - 1) + i - 1];
    }
  }
  for (int i = 1; i <= shared->dim_ICA_1; ++i) {
    for (int j = 1; j <= shared->dim_ICA_2; ++j) {
      shared->initial_ICA[i - 1 + shared->dim_ICA_1 * (j - 1)] = shared->init_ICA[shared->dim_init_ICA_1 * (j - 1) + i - 1];
    }
  }
  for (int i = 1; i <= shared->dim_ICM_1; ++i) {
    for (int j = 1; j <= shared->dim_ICM_2; ++j) {
      shared->initial_ICM[i - 1 + shared->dim_ICM_1 * (j - 1)] = shared->init_ICM[shared->dim_init_ICM_1 * (j - 1) + i - 1];
    }
  }
  for (int i = 1; i <= shared->dim_ID_1; ++i) {
    for (int j = 1; j <= shared->dim_ID_2; ++j) {
      shared->initial_ID[i - 1 + shared->dim_ID_1 * (j - 1)] = shared->init_ID[shared->dim_init_ID_1 * (j - 1) + i - 1];
    }
  }
  for (int i = 1; i <= shared->dim_P_1; ++i) {
    for (int j = 1; j <= shared->dim_P_2; ++j) {
      shared->initial_P[i - 1 + shared->dim_P_1 * (j - 1)] = shared->init_P[shared->dim_init_P_1 * (j - 1) + i - 1];
    }
  }
  for (int i = 1; i <= shared->dim_S_1; ++i) {
    for (int j = 1; j <= shared->dim_S_2; ++j) {
      shared->initial_S[i - 1 + shared->dim_S_1 * (j - 1)] = shared->init_S[shared->dim_init_S_1 * (j - 1) + i - 1];
    }
  }
  for (int i = 1; i <= shared->dim_T_1; ++i) {
    for (int j = 1; j <= shared->dim_T_2; ++j) {
      shared->initial_T[i - 1 + shared->dim_T_1 * (j - 1)] = shared->init_T[shared->dim_init_T_1 * (j - 1) + i - 1];
    }
  }
  for (int i = 1; i <= shared->dim_U_1; ++i) {
    for (int j = 1; j <= shared->dim_U_2; ++j) {
      shared->initial_U[i - 1 + shared->dim_U_1 * (j - 1)] = shared->init_U[shared->dim_init_U_1 * (j - 1) + i - 1];
    }
  }
  shared->offset_variable_A_preg = shared->dim_A + shared->dim_D + shared->dim_D_preg + shared->dim_FOIv + shared->dim_ince_delay + shared->dim_P + shared->dim_pregs + shared->dim_S + shared->dim_S_preg + shared->dim_T + shared->dim_T_preg + shared->dim_U + 4;
  shared->offset_variable_D_preg = shared->dim_A + shared->dim_D + shared->dim_FOIv + shared->dim_ince_delay + shared->dim_P + shared->dim_pregs + shared->dim_S + shared->dim_S_preg + shared->dim_T + shared->dim_T_preg + shared->dim_U + 4;
  shared->offset_variable_FOI = shared->dim_A + shared->dim_A_preg + shared->dim_D + shared->dim_D_preg + shared->dim_FOIv + shared->dim_IB + shared->dim_ICA + shared->dim_ICM + shared->dim_ID + shared->dim_ince_delay + shared->dim_P + shared->dim_P_preg + shared->dim_pregs + shared->dim_S + shared->dim_S_preg + shared->dim_T + shared->dim_T_preg + shared->dim_U + shared->dim_U_preg + 4;
  shared->offset_variable_IB = shared->dim_A + shared->dim_A_preg + shared->dim_D + shared->dim_D_preg + shared->dim_FOIv + shared->dim_ICA + shared->dim_ICM + shared->dim_ince_delay + shared->dim_P + shared->dim_P_preg + shared->dim_pregs + shared->dim_S + shared->dim_S_preg + shared->dim_T + shared->dim_T_preg + shared->dim_U + shared->dim_U_preg + 4;
  shared->offset_variable_ICA = shared->dim_A + shared->dim_A_preg + shared->dim_D + shared->dim_D_preg + shared->dim_FOIv + shared->dim_ICM + shared->dim_ince_delay + shared->dim_P + shared->dim_P_preg + shared->dim_pregs + shared->dim_S + shared->dim_S_preg + shared->dim_T + shared->dim_T_preg + shared->dim_U + shared->dim_U_preg + 4;
  shared->offset_variable_ICM = shared->dim_A + shared->dim_A_preg + shared->dim_D + shared->dim_D_preg + shared->dim_FOIv + shared->dim_ince_delay + shared->dim_P + shared->dim_P_preg + shared->dim_pregs + shared->dim_S + shared->dim_S_preg + shared->dim_T + shared->dim_T_preg + shared->dim_U + shared->dim_U_preg + 4;
  shared->offset_variable_ID = shared->dim_A + shared->dim_A_preg + shared->dim_D + shared->dim_D_preg + shared->dim_FOIv + shared->dim_IB + shared->dim_ICA + shared->dim_ICM + shared->dim_ince_delay + shared->dim_P + shared->dim_P_preg + shared->dim_pregs + shared->dim_S + shared->dim_S_preg + shared->dim_T + shared->dim_T_preg + shared->dim_U + shared->dim_U_preg + 4;
  shared->offset_variable_P_preg = shared->dim_A + shared->dim_A_preg + shared->dim_D + shared->dim_D_preg + shared->dim_FOIv + shared->dim_ince_delay + shared->dim_P + shared->dim_pregs + shared->dim_S + shared->dim_S_preg + shared->dim_T + shared->dim_T_preg + shared->dim_U + shared->dim_U_preg + 4;
  shared->offset_variable_T_preg = shared->dim_A + shared->dim_D + shared->dim_FOIv + shared->dim_ince_delay + shared->dim_P + shared->dim_pregs + shared->dim_S + shared->dim_S_preg + shared->dim_T + shared->dim_U + 4;
  shared->offset_variable_U_preg = shared->dim_A + shared->dim_A_preg + shared->dim_D + shared->dim_D_preg + shared->dim_FOIv + shared->dim_ince_delay + shared->dim_P + shared->dim_pregs + shared->dim_S + shared->dim_S_preg + shared->dim_T + shared->dim_T_preg + shared->dim_U + 4;
  for (int i = 1; i <= shared->preg_range; ++i) {
    for (int j = 1; j <= shared->nh; ++j) {
      shared->initial_A_preg[i - 1 + shared->dim_A_preg_1 * (j - 1)] = shared->init_A_preg[shared->dim_init_A_preg_1 * (j - 1) + i - 1];
    }
  }
  for (int i = 1; i <= shared->preg_range; ++i) {
    for (int j = 1; j <= shared->nh; ++j) {
      shared->initial_D_preg[i - 1 + shared->dim_D_preg_1 * (j - 1)] = shared->init_D_preg[shared->dim_init_D_preg_1 * (j - 1) + i - 1];
    }
  }
  for (int i = 1; i <= shared->dim_FOI_1; ++i) {
    for (int j = 1; j <= shared->dim_FOI_2; ++j) {
      for (int k = 1; k <= shared->dim_FOI_3; ++k) {
        shared->initial_FOI[i - 1 + shared->dim_FOI_1 * (j - 1) + shared->dim_FOI_12 * (k - 1)] = shared->init_FOI[shared->dim_init_FOI_12 * (k - 1) + shared->dim_init_FOI_1 * (j - 1) + i - 1];
      }
    }
  }
  for (int i = 1; i <= shared->preg_range; ++i) {
    for (int j = 1; j <= shared->nh; ++j) {
      shared->initial_P_preg[i - 1 + shared->dim_P_preg_1 * (j - 1)] = shared->init_P_preg[shared->dim_init_P_preg_1 * (j - 1) + i - 1];
    }
  }
  for (int i = 1; i <= shared->preg_range; ++i) {
    for (int j = 1; j <= shared->nh; ++j) {
      shared->initial_S_preg[i - 1 + shared->dim_S_preg_1 * (j - 1)] = shared->init_S_preg[shared->dim_init_S_preg_1 * (j - 1) + i - 1];
    }
  }
  for (int i = 1; i <= shared->preg_range; ++i) {
    for (int j = 1; j <= shared->nh; ++j) {
      shared->initial_T_preg[i - 1 + shared->dim_T_preg_1 * (j - 1)] = shared->init_T_preg[shared->dim_init_T_preg_1 * (j - 1) + i - 1];
    }
  }
  for (int i = 1; i <= shared->preg_range; ++i) {
    for (int j = 1; j <= shared->nh; ++j) {
      shared->initial_U_preg[i - 1 + shared->dim_U_preg_1 * (j - 1)] = shared->init_U_preg[shared->dim_init_U_preg_1 * (j - 1) + i - 1];
    }
  }
  return mode::pars_type<mipodinmodel>(shared, internal);
}
template <>
cpp11::sexp mode_info<mipodinmodel>(const mode::pars_type<mipodinmodel>& pars) {
  const mipodinmodel::internal_type internal = pars.internal;
  const std::shared_ptr<const mipodinmodel::shared_type> shared = pars.shared;
  cpp11::writable::strings nms({"betaa_td", "Sv", "Ev", "Iv", "pregs", "FOIv", "ince_delay", "S", "T", "D", "A", "U", "P", "S_preg", "T_preg", "D_preg", "A_preg", "U_preg", "P_preg", "ICM", "ICA", "IB", "ID", "FOI", "ICM_pre", "Sout", "Tout", "Dout", "Aout", "Uout", "Pout", "S_preg_out", "T_preg_out", "D_preg_out", "A_preg_out", "U_preg_out", "P_preg_out", "Y_preg_out", "prev", "inc05", "inc", "prev_rdt_primi"});
  cpp11::writable::list dim(42);
  dim[0] = cpp11::writable::integers({1});
  dim[1] = cpp11::writable::integers({1});
  dim[2] = cpp11::writable::integers({1});
  dim[3] = cpp11::writable::integers({1});
  dim[4] = cpp11::writable::integers({shared->dim_pregs});
  dim[5] = cpp11::writable::integers({shared->dim_FOIv});
  dim[6] = cpp11::writable::integers({shared->dim_ince_delay});
  dim[7] = cpp11::writable::integers({shared->dim_S_1, shared->dim_S_2});
  dim[8] = cpp11::writable::integers({shared->dim_T_1, shared->dim_T_2});
  dim[9] = cpp11::writable::integers({shared->dim_D_1, shared->dim_D_2});
  dim[10] = cpp11::writable::integers({shared->dim_A_1, shared->dim_A_2});
  dim[11] = cpp11::writable::integers({shared->dim_U_1, shared->dim_U_2});
  dim[12] = cpp11::writable::integers({shared->dim_P_1, shared->dim_P_2});
  dim[13] = cpp11::writable::integers({shared->dim_S_preg_1, shared->dim_S_preg_2});
  dim[14] = cpp11::writable::integers({shared->dim_T_preg_1, shared->dim_T_preg_2});
  dim[15] = cpp11::writable::integers({shared->dim_D_preg_1, shared->dim_D_preg_2});
  dim[16] = cpp11::writable::integers({shared->dim_A_preg_1, shared->dim_A_preg_2});
  dim[17] = cpp11::writable::integers({shared->dim_U_preg_1, shared->dim_U_preg_2});
  dim[18] = cpp11::writable::integers({shared->dim_P_preg_1, shared->dim_P_preg_2});
  dim[19] = cpp11::writable::integers({shared->dim_ICM_1, shared->dim_ICM_2});
  dim[20] = cpp11::writable::integers({shared->dim_ICA_1, shared->dim_ICA_2});
  dim[21] = cpp11::writable::integers({shared->dim_IB_1, shared->dim_IB_2});
  dim[22] = cpp11::writable::integers({shared->dim_ID_1, shared->dim_ID_2});
  dim[23] = cpp11::writable::integers({shared->dim_FOI_1, shared->dim_FOI_2, shared->dim_FOI_3});
  dim[24] = cpp11::writable::integers({1});
  dim[25] = cpp11::writable::integers({1});
  dim[26] = cpp11::writable::integers({1});
  dim[27] = cpp11::writable::integers({1});
  dim[28] = cpp11::writable::integers({1});
  dim[29] = cpp11::writable::integers({1});
  dim[30] = cpp11::writable::integers({1});
  dim[31] = cpp11::writable::integers({1});
  dim[32] = cpp11::writable::integers({1});
  dim[33] = cpp11::writable::integers({1});
  dim[34] = cpp11::writable::integers({1});
  dim[35] = cpp11::writable::integers({1});
  dim[36] = cpp11::writable::integers({1});
  dim[37] = cpp11::writable::integers({1});
  dim[38] = cpp11::writable::integers({1});
  dim[39] = cpp11::writable::integers({1});
  dim[40] = cpp11::writable::integers({1});
  dim[41] = cpp11::writable::integers({1});
  dim.names() = nms;
  cpp11::writable::list index(42);
  index[0] = cpp11::writable::integers({1});
  index[1] = cpp11::writable::integers({2});
  index[2] = cpp11::writable::integers({3});
  index[3] = cpp11::writable::integers({4});
  index[4] = integer_sequence(5, shared->dim_pregs);
  index[5] = integer_sequence(shared->offset_variable_FOIv + 1, shared->dim_FOIv);
  index[6] = integer_sequence(shared->offset_variable_ince_delay + 1, shared->dim_ince_delay);
  index[7] = integer_sequence(shared->offset_variable_S + 1, shared->dim_S);
  index[8] = integer_sequence(shared->offset_variable_T + 1, shared->dim_T);
  index[9] = integer_sequence(shared->offset_variable_D + 1, shared->dim_D);
  index[10] = integer_sequence(shared->offset_variable_A + 1, shared->dim_A);
  index[11] = integer_sequence(shared->offset_variable_U + 1, shared->dim_U);
  index[12] = integer_sequence(shared->offset_variable_P + 1, shared->dim_P);
  index[13] = integer_sequence(shared->offset_variable_S_preg + 1, shared->dim_S_preg);
  index[14] = integer_sequence(shared->offset_variable_T_preg + 1, shared->dim_T_preg);
  index[15] = integer_sequence(shared->offset_variable_D_preg + 1, shared->dim_D_preg);
  index[16] = integer_sequence(shared->offset_variable_A_preg + 1, shared->dim_A_preg);
  index[17] = integer_sequence(shared->offset_variable_U_preg + 1, shared->dim_U_preg);
  index[18] = integer_sequence(shared->offset_variable_P_preg + 1, shared->dim_P_preg);
  index[19] = integer_sequence(shared->offset_variable_ICM + 1, shared->dim_ICM);
  index[20] = integer_sequence(shared->offset_variable_ICA + 1, shared->dim_ICA);
  index[21] = integer_sequence(shared->offset_variable_IB + 1, shared->dim_IB);
  index[22] = integer_sequence(shared->offset_variable_ID + 1, shared->dim_ID);
  index[23] = integer_sequence(shared->offset_variable_FOI + 1, shared->dim_FOI);
  index[24] = cpp11::writable::integers({0 + shared->offset_variable_FOI + 1 + shared->dim_FOI});
  index[25] = cpp11::writable::integers({1 + shared->offset_variable_FOI + 1 + shared->dim_FOI});
  index[26] = cpp11::writable::integers({2 + shared->offset_variable_FOI + 1 + shared->dim_FOI});
  index[27] = cpp11::writable::integers({3 + shared->offset_variable_FOI + 1 + shared->dim_FOI});
  index[28] = cpp11::writable::integers({4 + shared->offset_variable_FOI + 1 + shared->dim_FOI});
  index[29] = cpp11::writable::integers({5 + shared->offset_variable_FOI + 1 + shared->dim_FOI});
  index[30] = cpp11::writable::integers({6 + shared->offset_variable_FOI + 1 + shared->dim_FOI});
  index[31] = cpp11::writable::integers({7 + shared->offset_variable_FOI + 1 + shared->dim_FOI});
  index[32] = cpp11::writable::integers({8 + shared->offset_variable_FOI + 1 + shared->dim_FOI});
  index[33] = cpp11::writable::integers({9 + shared->offset_variable_FOI + 1 + shared->dim_FOI});
  index[34] = cpp11::writable::integers({10 + shared->offset_variable_FOI + 1 + shared->dim_FOI});
  index[35] = cpp11::writable::integers({11 + shared->offset_variable_FOI + 1 + shared->dim_FOI});
  index[36] = cpp11::writable::integers({12 + shared->offset_variable_FOI + 1 + shared->dim_FOI});
  index[37] = cpp11::writable::integers({13 + shared->offset_variable_FOI + 1 + shared->dim_FOI});
  index[38] = cpp11::writable::integers({14 + shared->offset_variable_FOI + 1 + shared->dim_FOI});
  index[39] = cpp11::writable::integers({15 + shared->offset_variable_FOI + 1 + shared->dim_FOI});
  index[40] = cpp11::writable::integers({16 + shared->offset_variable_FOI + 1 + shared->dim_FOI});
  index[41] = cpp11::writable::integers({17 + shared->offset_variable_FOI + 1 + shared->dim_FOI});
  index.names() = nms;
  size_t len = shared->offset_variable_FOI + shared->dim_FOI + 18;
  using namespace cpp11::literals;
  return cpp11::writable::list({
           "dim"_nm = dim,
           "len"_nm = len,
           "index"_nm = index});
}
}

[[cpp11::register]]
SEXP mode_mipodinmodel_alloc(cpp11::list r_pars, double time,
                         size_t n_particles, size_t n_threads,
                         cpp11::sexp control, cpp11::sexp seed) {
  return mode::r::mode_alloc<mipodinmodel>(r_pars, time, n_particles, n_threads,
      control, seed);
}

[[cpp11::register]]
cpp11::sexp mode_mipodinmodel_control(SEXP ptr) {
  return mode::r::mode_control<mode::container<mipodinmodel>>(ptr);
}

[[cpp11::register]]
double mode_mipodinmodel_time(SEXP ptr) {
  return mode::r::mode_time<mode::container<mipodinmodel>>(ptr);
}

[[cpp11::register]]
cpp11::sexp mode_mipodinmodel_run(SEXP ptr, double end_time) {
  return mode::r::mode_run<mode::container<mipodinmodel>>(ptr, end_time);
}

[[cpp11::register]]
cpp11::sexp mode_mipodinmodel_state_full(SEXP ptr) {
  return mode::r::mode_state_full<mode::container<mipodinmodel>>(ptr);
}

[[cpp11::register]]
cpp11::sexp mode_mipodinmodel_state(SEXP ptr, SEXP index) {
  return mode::r::mode_state<mode::container<mipodinmodel>>(ptr, index);
}

[[cpp11::register]]
cpp11::sexp mode_mipodinmodel_stats(SEXP ptr) {
  return mode::r::mode_stats<mode::container<mipodinmodel>>(ptr);
}

[[cpp11::register]]
cpp11::sexp mode_mipodinmodel_update_state(SEXP ptr,
                                       SEXP pars,
                                       SEXP time,
                                       SEXP state,
                                       SEXP index,
                                       SEXP set_initial_state,
                                       SEXP reset_step_size) {
  return mode::r::mode_update_state<mipodinmodel>(ptr,
      pars, time, state, index, set_initial_state, reset_step_size);
}

[[cpp11::register]]
void mode_mipodinmodel_set_index(SEXP ptr, SEXP index) {
  return mode::r::mode_set_index<mode::container<mipodinmodel>>(ptr, index);
}

[[cpp11::register]]
void mode_mipodinmodel_set_stochastic_schedule(SEXP ptr, SEXP time) {
  return mode::r::mode_set_stochastic_schedule<mode::container<mipodinmodel>>(ptr, time);
}

[[cpp11::register]]
size_t mode_mipodinmodel_n_variables(SEXP ptr) {
  return mode::r::mode_n_variables<mode::container<mipodinmodel>>(ptr);
}

[[cpp11::register]]
size_t mode_mipodinmodel_n_state_run(SEXP ptr) {
  return mode::r::mode_n_state_run<mode::container<mipodinmodel>>(ptr);
}

[[cpp11::register]]
size_t mode_mipodinmodel_n_state_full(SEXP ptr) {
  return mode::r::mode_n_state_full<mode::container<mipodinmodel>>(ptr);
}

[[cpp11::register]]
void mode_mipodinmodel_reorder(SEXP ptr, SEXP index) {
  return mode::r::mode_reorder<mode::container<mipodinmodel>>(ptr, index);
}

[[cpp11::register]]
void mode_mipodinmodel_set_n_threads(SEXP ptr, int n_threads) {
  return mode::r::mode_set_n_threads<mode::container<mipodinmodel>>(ptr, n_threads);
}

[[cpp11::register]]
bool mode_mipodinmodel_has_openmp() {
#ifdef _OPENMP
  return true;
#else
  return false;
#endif
}
