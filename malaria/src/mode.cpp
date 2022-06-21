#include <mode/r/mode.hpp>

template <typename real_type, typename container>
__host__ __device__ real_type odin_sum1(const container x, size_t from, size_t to);
// [[dust::class(odinmodel)]]
// [[dust::param(init_Ih, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(init_Iv, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(init_Sv, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(nrates, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(a, has_default = TRUE, default_value = list("/", 1L, 3L), rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(r, has_default = TRUE, default_value = list("/", 1L, 100L), rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(tau, has_default = TRUE, default_value = 12L, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
class odinmodel {
public:
  using real_type = double;
  using rng_state_type = dust::random::generator<real_type>;
  using data_type = mode::no_data;
  struct shared_type {
    real_type a;
    real_type bh;
    real_type bv;
    int dim_Ev;
    real_type init_Ih;
    real_type init_Iv;
    real_type init_Sh;
    real_type init_Sv;
    real_type initial_beta;
    std::vector<real_type> initial_Ev;
    real_type initial_Ih;
    real_type initial_Iv;
    real_type initial_Sh;
    real_type initial_Sv;
    real_type mu;
    int nrates;
    real_type p;
    real_type r;
    real_type tau;
  };
  struct internal_type {
  };
  odinmodel(const mode::pars_type<odinmodel>& pars) :
    shared(pars.shared), internal(pars.internal) {
  }
  size_t n_variables() {
    return shared->dim_Ev + 5;
  }
  size_t n_output() {
    return 3;
  }
  std::vector<real_type> initial(size_t t) {
    std::vector<real_type> state(shared->dim_Ev + 5);
    state[0] = shared->initial_Sh;
    state[1] = shared->initial_Ih;
    state[2] = shared->initial_Sv;
    state[3] = shared->initial_Iv;
    state[4] = shared->initial_beta;
    std::copy(shared->initial_Ev.begin(), shared->initial_Ev.end(), state.begin() + 5);
    return state;
  }
  void update_stochastic(double t,
                         const std::vector<double>& state,
                         rng_state_type& rng_state,
                         std::vector<double>& state_next) {
    const int beta_index = 4;
    const double beta_volatility = 0.5;
    state_next[beta_index] *= std::exp(dust::random::normal<double>(rng_state, 0, beta_volatility));
  }
  void rhs(double t, const std::vector<double>& state, std::vector<double>& dstatedt) {
    const real_type Sh = state[0];
    const real_type Ih = state[1];
    const real_type Sv = state[2];
    const real_type * Ev = state.data() + 5;
    const real_type Iv = state[3];
    const real_type beta = state[4];
    dstatedt[4] = 0;
    real_type foi_h = shared->a * shared->bh * Iv;
    real_type foi_v = shared->a * shared->bv * Ih;
    {
       int i = 1;
       dstatedt[5 + i - 1] = foi_v * Sv - (shared->nrates / (real_type) shared->tau) * Ev[i - 1] - shared->mu * Ev[i - 1];
    }
    for (int i = 2; i <= shared->nrates; ++i) {
      dstatedt[5 + i - 1] = (shared->nrates / (real_type) shared->tau) * Ev[i - 1 - 1] - (shared->nrates / (real_type) shared->tau) * Ev[i - 1] - shared->mu * Ev[i - 1];
    }
    dstatedt[1] = foi_h * Sh - shared->r * Ih;
    dstatedt[3] = (shared->nrates / (real_type) shared->tau) * Ev[shared->nrates - 1] - shared->mu * Iv;
    dstatedt[0] = - foi_h * Sh + shared->r * Ih;
    dstatedt[2] = beta * shared->init_Sv - foi_v * Sv - shared->mu * Sv;
  }
  void output(double t, const std::vector<double>& state, std::vector<double>& output) {
    const real_type Sh = state[0];
    const real_type Ih = state[1];
    const real_type Sv = state[2];
    const real_type * Ev = state.data() + 5;
    const real_type Iv = state[3];
    real_type N = Sh + Ih;
    output[2] = odin_sum1<real_type>(Ev, 0, shared->dim_Ev);
    output[0] = Ih / (real_type) N;
    real_type V = Sv + odin_sum1<real_type>(Ev, 0, shared->dim_Ev) + Iv;
    output[1] = Iv / (real_type) V;
  }
private:
  std::shared_ptr<const shared_type> shared;
  internal_type internal;
};
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
mode::pars_type<odinmodel> mode_pars<odinmodel>(cpp11::list user) {
  using real_type = typename odinmodel::real_type;
  auto shared = std::make_shared<odinmodel::shared_type>();
  odinmodel::internal_type internal;
  shared->bh = 0.050000000000000003;
  shared->bv = 0.050000000000000003;
  shared->p = 0.90000000000000002;
  shared->initial_beta = - std::log(shared->p);
  shared->mu = - std::log(shared->p);
  shared->init_Ih = NA_REAL;
  shared->init_Iv = NA_REAL;
  shared->init_Sv = NA_REAL;
  shared->nrates = NA_INTEGER;
  shared->a = 1 / (real_type) 3;
  shared->r = 1 / (real_type) 100;
  shared->tau = 12;
  shared->a = user_get_scalar<real_type>(user, "a", shared->a, NA_REAL, NA_REAL);
  shared->init_Ih = user_get_scalar<real_type>(user, "init_Ih", shared->init_Ih, NA_REAL, NA_REAL);
  shared->init_Iv = user_get_scalar<real_type>(user, "init_Iv", shared->init_Iv, NA_REAL, NA_REAL);
  shared->init_Sv = user_get_scalar<real_type>(user, "init_Sv", shared->init_Sv, NA_REAL, NA_REAL);
  shared->nrates = user_get_scalar<int>(user, "nrates", shared->nrates, NA_INTEGER, NA_INTEGER);
  shared->r = user_get_scalar<real_type>(user, "r", shared->r, NA_REAL, NA_REAL);
  shared->tau = user_get_scalar<real_type>(user, "tau", shared->tau, NA_REAL, NA_REAL);
  shared->dim_Ev = shared->nrates;
  shared->init_Sh = 1 - shared->init_Ih;
  shared->initial_Ih = shared->init_Ih;
  shared->initial_Iv = shared->init_Iv;
  shared->initial_Sv = shared->init_Sv;
  shared->initial_Ev = std::vector<real_type>(shared->dim_Ev);
  for (int i = 1; i <= shared->dim_Ev; ++i) {
    shared->initial_Ev[i - 1] = 0;
  }
  shared->initial_Sh = shared->init_Sh;
  return mode::pars_type<odinmodel>(shared, internal);
}
template <>
cpp11::sexp mode_info<odinmodel>(const mode::pars_type<odinmodel>& pars) {
  const odinmodel::internal_type internal = pars.internal;
  const std::shared_ptr<const odinmodel::shared_type> shared = pars.shared;
  cpp11::writable::strings nms({"Sh", "Ih", "Sv", "Iv", "beta", "Ev"});
  cpp11::writable::list dim(6);
  dim[0] = cpp11::writable::integers({1});
  dim[1] = cpp11::writable::integers({1});
  dim[2] = cpp11::writable::integers({1});
  dim[3] = cpp11::writable::integers({1});
  dim[4] = cpp11::writable::integers({1});
  dim[5] = cpp11::writable::integers({shared->dim_Ev});
  dim.names() = nms;
  cpp11::writable::list index(6);
  index[0] = cpp11::writable::integers({1});
  index[1] = cpp11::writable::integers({2});
  index[2] = cpp11::writable::integers({3});
  index[3] = cpp11::writable::integers({4});
  index[4] = cpp11::writable::integers({5});
  index[5] = integer_sequence(6, shared->dim_Ev);
  index.names() = nms;
  size_t len = 5 + shared->dim_Ev;
  using namespace cpp11::literals;
  return cpp11::writable::list({
           "dim"_nm = dim,
           "len"_nm = len,
           "index"_nm = index});
}
}

[[cpp11::register]]
SEXP mode_odinmodel_alloc(cpp11::list r_pars, double time,
                         size_t n_particles, size_t n_threads,
                         cpp11::sexp control, cpp11::sexp seed) {
  return mode::r::mode_alloc<odinmodel>(r_pars, time, n_particles, n_threads,
      control, seed);
}

[[cpp11::register]]
cpp11::sexp mode_odinmodel_control(SEXP ptr) {
  return mode::r::mode_control<mode::container<odinmodel>>(ptr);
}

[[cpp11::register]]
double mode_odinmodel_time(SEXP ptr) {
  return mode::r::mode_time<mode::container<odinmodel>>(ptr);
}

[[cpp11::register]]
cpp11::sexp mode_odinmodel_run(SEXP ptr, double end_time) {
  return mode::r::mode_run<mode::container<odinmodel>>(ptr, end_time);
}

[[cpp11::register]]
cpp11::sexp mode_odinmodel_state_full(SEXP ptr) {
  return mode::r::mode_state_full<mode::container<odinmodel>>(ptr);
}

[[cpp11::register]]
cpp11::sexp mode_odinmodel_state(SEXP ptr, SEXP index) {
  return mode::r::mode_state<mode::container<odinmodel>>(ptr, index);
}

[[cpp11::register]]
cpp11::sexp mode_odinmodel_stats(SEXP ptr) {
  return mode::r::mode_stats<mode::container<odinmodel>>(ptr);
}

[[cpp11::register]]
cpp11::sexp mode_odinmodel_update_state(SEXP ptr,
                                       SEXP pars,
                                       SEXP time,
                                       SEXP state,
                                       SEXP index,
                                       SEXP set_initial_state,
                                       SEXP reset_step_size) {
  return mode::r::mode_update_state<odinmodel>(ptr,
      pars, time, state, index, set_initial_state, reset_step_size);
}

[[cpp11::register]]
void mode_odinmodel_set_index(SEXP ptr, SEXP index) {
  return mode::r::mode_set_index<mode::container<odinmodel>>(ptr, index);
}

[[cpp11::register]]
void mode_odinmodel_set_stochastic_schedule(SEXP ptr, SEXP time) {
  return mode::r::mode_set_stochastic_schedule<mode::container<odinmodel>>(ptr, time);
}

[[cpp11::register]]
size_t mode_odinmodel_n_variables(SEXP ptr) {
  return mode::r::mode_n_variables<mode::container<odinmodel>>(ptr);
}

[[cpp11::register]]
size_t mode_odinmodel_n_state_run(SEXP ptr) {
  return mode::r::mode_n_state_run<mode::container<odinmodel>>(ptr);
}

[[cpp11::register]]
size_t mode_odinmodel_n_state_full(SEXP ptr) {
  return mode::r::mode_n_state_full<mode::container<odinmodel>>(ptr);
}

[[cpp11::register]]
void mode_odinmodel_reorder(SEXP ptr, SEXP index) {
  return mode::r::mode_reorder<mode::container<odinmodel>>(ptr, index);
}

[[cpp11::register]]
void mode_odinmodel_set_n_threads(SEXP ptr, int n_threads) {
  return mode::r::mode_set_n_threads<mode::container<odinmodel>>(ptr, n_threads);
}

[[cpp11::register]]
bool mode_odinmodel_has_openmp() {
#ifdef _OPENMP
  return true;
#else
  return false;
#endif
}
