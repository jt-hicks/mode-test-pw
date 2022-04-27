class malaria {
public:
  using data_type = mode::no_data;
  using internal_type = mode::no_internal;
  using rng_state_type = dust::random::generator<double>;

  struct shared_type {
    int nrates;
    double a;
    double r;
    double tau;
    double bh;
    double bv;
    double mu;
    double init_Ih;
    double init_Sv;
    double init_Iv;
    double init_beta;
    double beta_volatility;
  };

  malaria(const mode::pars_type<malaria>& pars) : shared(pars.shared) {
  }

  void rhs(double t,
           const std::vector<double>& y,
           std::vector<double>& dydt) {
    const double Sh = y[0], Ih = y[1], Sv = y[2], Iv = y[3], beta = y[4];
    const double* Ev = y.data() + 5;
    double * deriv_Ev = dydt.data() + 5;

    const double progress_Ev = shared->nrates / shared->tau;
    const double foi_h = shared->a * shared->bh * Iv;
    const double foi_v = shared->a * shared->bv * Ih;

    dydt[0] = -foi_h * Sh + shared->r * Ih; // Sh
    dydt[1] =  foi_h * Sh - shared->r * Ih; // Ih
    dydt[2] = beta * shared->init_Sv - foi_v * Sv - shared->mu * Sv; // Sv
    dydt[3] = progress_Ev * Ev[shared->nrates - 1] - shared->mu * Iv;
    dydt[4] = 0; // beta

    deriv_Ev[0] = foi_v * Sv - progress_Ev * Ev[0] - shared->mu * Ev[0];
    for (int i = 1; i < shared->nrates; ++i) {
      deriv_Ev[i] = progress_Ev * Ev[i - 1] -
        progress_Ev * Ev[i] - shared->mu * Ev[i];
    }
  }

  void update_stochastic(double t, std::vector<double>& y,
                         rng_state_type& rng_state) {
    const double r = dust::random::normal<double>(rng_state, 0,
                                                  shared->beta_volatility);
    y[4] *= std::exp(r); // updates beta
  }

  std::vector<double> initial(double time) {
    // ordering: Sh, Ih, Sv, Iv, beta, Ev[]
    std::vector<double> ret(size(), 0);
    ret[0] = 1 - shared->init_Ih; // Sh
    ret[1] = shared->init_Ih;     // Ih
    ret[2] = shared->init_Sv;     // Sv
    ret[3] = shared->init_Iv;     // Ev
    ret[4] = shared->init_beta;   // beta
    return ret;
  }

  size_t size() const {
    return 5 + shared->nrates;
  }

private:
  mode::shared_ptr<malaria> shared;
};

template <typename T>
inline T with_default(cpp11::sexp value, T default_value) {
  return value == R_NilValue ? default_value : cpp11::as_cpp<T>(value);
}

namespace mode {

template <>
mode::pars_type<malaria> mode_pars<malaria>(cpp11::list pars) {
  int nrates = with_default<int>(pars["nrates"], 15);
  double a = with_default<double>(pars["a"], 1.0 / 3.0);
  double r = with_default<double>(pars["r"], 1.0 / 100.0);
  double tau = with_default<double>(pars["tau"], 12);
  double init_Ih = with_default<double>(pars["init_Ih"], 0.8);
  double init_Sv = with_default<double>(pars["init_Sv"], 100);
  double init_Iv = with_default<double>(pars["init_Iv"], 1);
  double beta_volatility = with_default<double>(pars["beta_volatility"], 0.5);
  double bh = 0.05;
  double bv = 0.05;
  double mu = -std::log(0.9);
  double init_beta = mu;

  malaria::shared_type shared{nrates, a, r, tau, bh, bv, mu,
                              init_Ih, init_Sv, init_Iv,
                              init_beta, beta_volatility};
  return mode::pars_type<malaria>(shared);
}

}
