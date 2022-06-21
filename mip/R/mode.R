mipodinmodel <- R6::R6Class(
  "mode",
  cloneable = FALSE,

  private = list(
    pars_ = NULL,
    info_ = NULL,
    ptr_ = NULL,
    n_particles_ = NULL,
    n_threads_ = NULL,
    index_ = NULL
  ),

  public = list(
    initialize = function(pars, time, n_particles, n_threads = 1L,
                          control = NULL, seed = NULL) {
      res <- mode_mipodinmodel_alloc(pars, time, n_particles,
                                 n_threads, control, seed)
      private$ptr_ <- res[[1]]
      private$info_ <- res[[2]]
      private$pars_ <- pars
      private$n_particles_ <- n_particles
      private$n_threads_ <- n_threads
    },

    name = function() {
      "mipodinmodel"
    },

    n_particles = function() {
       private$n_particles_
    },

    pars = function() {
      private$pars_
    },

    info = function() {
      private$info_
    },

    time = function() {
      mode_mipodinmodel_time(private$ptr_)
    },

    set_index = function(index) {
      mode_mipodinmodel_set_index(private$ptr_, index)
      private$index_ <- index
      invisible()
    },

    set_stochastic_schedule = function(time) {
      mode_mipodinmodel_set_stochastic_schedule(private$ptr_, time)
      invisible()
    },

    index = function() {
      private$index_
    },

    control = function() {
      mode_mipodinmodel_control(private$ptr_)
    },

    run = function(end_time) {
      m <- mode_mipodinmodel_run(private$ptr_, end_time)
      rownames(m) <- names(private$index_)
      m
    },

    statistics = function() {
      mode_mipodinmodel_stats(private$ptr_)
    },

    n_state_run = function() {
      mode_mipodinmodel_n_state_run(private$ptr_)
    },

    n_state_full = function() {
      mode_mipodinmodel_n_state_full(private$ptr_)
    },

    n_threads = function() {
      private$n_threads_
    },

    set_n_threads = function(n_threads) {
      prev <- private$n_threads_
      mode_mipodinmodel_set_n_threads(private$ptr_, n_threads)
      private$n_threads_ <- n_threads
      invisible(prev)
    },

    update_state = function(pars = NULL, time = NULL,
                            state = NULL, index = NULL,
                            set_initial_state = NULL, reset_step_size = NULL) {
      info <- mode_mipodinmodel_update_state(private$ptr_, pars, time, state, index,
                                 set_initial_state, reset_step_size)
      if (!is.null(pars)) {
        private$pars_ <- pars
        private$info_ <- info
      }
    },

    reorder = function(index) {
      mode_mipodinmodel_reorder(private$ptr_, index)
      invisible()
    },

    state = function(index = NULL) {
      if (is.null(index)) {
        mode_mipodinmodel_state_full(private$ptr_)
      } else {
        mode_mipodinmodel_state(private$ptr_, index)
      }
    },

    has_openmp = function() {
      mode_mipodinmodel_has_openmp()
    }

  ))
