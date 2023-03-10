# Generated by cpp11: do not edit by hand

dust_seir_capabilities <- function() {
  .Call(`_seir_dust_seir_capabilities`)
}

dust_seir_gpu_info <- function() {
  .Call(`_seir_dust_seir_gpu_info`)
}

dust_cpu_seir_alloc <- function(r_pars, pars_multi, step, r_n_particles, n_threads, r_seed, deterministic, gpu_config) {
  .Call(`_seir_dust_cpu_seir_alloc`, r_pars, pars_multi, step, r_n_particles, n_threads, r_seed, deterministic, gpu_config)
}

dust_cpu_seir_run <- function(ptr, step_end) {
  .Call(`_seir_dust_cpu_seir_run`, ptr, step_end)
}

dust_cpu_seir_simulate <- function(ptr, step_end) {
  .Call(`_seir_dust_cpu_seir_simulate`, ptr, step_end)
}

dust_cpu_seir_set_index <- function(ptr, r_index) {
  .Call(`_seir_dust_cpu_seir_set_index`, ptr, r_index)
}

dust_cpu_seir_update_state <- function(ptr, r_pars, r_state, r_step, r_set_initial_state) {
  .Call(`_seir_dust_cpu_seir_update_state`, ptr, r_pars, r_state, r_step, r_set_initial_state)
}

dust_cpu_seir_state <- function(ptr, r_index) {
  .Call(`_seir_dust_cpu_seir_state`, ptr, r_index)
}

dust_cpu_seir_step <- function(ptr) {
  .Call(`_seir_dust_cpu_seir_step`, ptr)
}

dust_cpu_seir_reorder <- function(ptr, r_index) {
  invisible(.Call(`_seir_dust_cpu_seir_reorder`, ptr, r_index))
}

dust_cpu_seir_resample <- function(ptr, r_weights) {
  .Call(`_seir_dust_cpu_seir_resample`, ptr, r_weights)
}

dust_cpu_seir_rng_state <- function(ptr, first_only, last_only) {
  .Call(`_seir_dust_cpu_seir_rng_state`, ptr, first_only, last_only)
}

dust_cpu_seir_set_rng_state <- function(ptr, rng_state) {
  .Call(`_seir_dust_cpu_seir_set_rng_state`, ptr, rng_state)
}

dust_cpu_seir_set_data <- function(ptr, data) {
  .Call(`_seir_dust_cpu_seir_set_data`, ptr, data)
}

dust_cpu_seir_compare_data <- function(ptr) {
  .Call(`_seir_dust_cpu_seir_compare_data`, ptr)
}

dust_cpu_seir_filter <- function(ptr, step_end, save_trajectories, step_snapshot, min_log_likelihood) {
  .Call(`_seir_dust_cpu_seir_filter`, ptr, step_end, save_trajectories, step_snapshot, min_log_likelihood)
}

dust_cpu_seir_set_n_threads <- function(ptr, n_threads) {
  invisible(.Call(`_seir_dust_cpu_seir_set_n_threads`, ptr, n_threads))
}

dust_cpu_seir_n_state <- function(ptr) {
  .Call(`_seir_dust_cpu_seir_n_state`, ptr)
}

dust_seirdage_capabilities <- function() {
  .Call(`_seir_dust_seirdage_capabilities`)
}

dust_seirdage_gpu_info <- function() {
  .Call(`_seir_dust_seirdage_gpu_info`)
}

dust_cpu_seirdage_alloc <- function(r_pars, pars_multi, step, r_n_particles, n_threads, r_seed, deterministic, gpu_config) {
  .Call(`_seir_dust_cpu_seirdage_alloc`, r_pars, pars_multi, step, r_n_particles, n_threads, r_seed, deterministic, gpu_config)
}

dust_cpu_seirdage_run <- function(ptr, step_end) {
  .Call(`_seir_dust_cpu_seirdage_run`, ptr, step_end)
}

dust_cpu_seirdage_simulate <- function(ptr, step_end) {
  .Call(`_seir_dust_cpu_seirdage_simulate`, ptr, step_end)
}

dust_cpu_seirdage_set_index <- function(ptr, r_index) {
  .Call(`_seir_dust_cpu_seirdage_set_index`, ptr, r_index)
}

dust_cpu_seirdage_update_state <- function(ptr, r_pars, r_state, r_step, r_set_initial_state) {
  .Call(`_seir_dust_cpu_seirdage_update_state`, ptr, r_pars, r_state, r_step, r_set_initial_state)
}

dust_cpu_seirdage_state <- function(ptr, r_index) {
  .Call(`_seir_dust_cpu_seirdage_state`, ptr, r_index)
}

dust_cpu_seirdage_step <- function(ptr) {
  .Call(`_seir_dust_cpu_seirdage_step`, ptr)
}

dust_cpu_seirdage_reorder <- function(ptr, r_index) {
  invisible(.Call(`_seir_dust_cpu_seirdage_reorder`, ptr, r_index))
}

dust_cpu_seirdage_resample <- function(ptr, r_weights) {
  .Call(`_seir_dust_cpu_seirdage_resample`, ptr, r_weights)
}

dust_cpu_seirdage_rng_state <- function(ptr, first_only, last_only) {
  .Call(`_seir_dust_cpu_seirdage_rng_state`, ptr, first_only, last_only)
}

dust_cpu_seirdage_set_rng_state <- function(ptr, rng_state) {
  .Call(`_seir_dust_cpu_seirdage_set_rng_state`, ptr, rng_state)
}

dust_cpu_seirdage_set_data <- function(ptr, data) {
  .Call(`_seir_dust_cpu_seirdage_set_data`, ptr, data)
}

dust_cpu_seirdage_compare_data <- function(ptr) {
  .Call(`_seir_dust_cpu_seirdage_compare_data`, ptr)
}

dust_cpu_seirdage_filter <- function(ptr, step_end, save_trajectories, step_snapshot, min_log_likelihood) {
  .Call(`_seir_dust_cpu_seirdage_filter`, ptr, step_end, save_trajectories, step_snapshot, min_log_likelihood)
}

dust_cpu_seirdage_set_n_threads <- function(ptr, n_threads) {
  invisible(.Call(`_seir_dust_cpu_seirdage_set_n_threads`, ptr, n_threads))
}

dust_cpu_seirdage_n_state <- function(ptr) {
  .Call(`_seir_dust_cpu_seirdage_n_state`, ptr)
}

dust_seirhdage_capabilities <- function() {
  .Call(`_seir_dust_seirhdage_capabilities`)
}

dust_seirhdage_gpu_info <- function() {
  .Call(`_seir_dust_seirhdage_gpu_info`)
}

dust_cpu_seirhdage_alloc <- function(r_pars, pars_multi, step, r_n_particles, n_threads, r_seed, deterministic, gpu_config) {
  .Call(`_seir_dust_cpu_seirhdage_alloc`, r_pars, pars_multi, step, r_n_particles, n_threads, r_seed, deterministic, gpu_config)
}

dust_cpu_seirhdage_run <- function(ptr, step_end) {
  .Call(`_seir_dust_cpu_seirhdage_run`, ptr, step_end)
}

dust_cpu_seirhdage_simulate <- function(ptr, step_end) {
  .Call(`_seir_dust_cpu_seirhdage_simulate`, ptr, step_end)
}

dust_cpu_seirhdage_set_index <- function(ptr, r_index) {
  .Call(`_seir_dust_cpu_seirhdage_set_index`, ptr, r_index)
}

dust_cpu_seirhdage_update_state <- function(ptr, r_pars, r_state, r_step, r_set_initial_state) {
  .Call(`_seir_dust_cpu_seirhdage_update_state`, ptr, r_pars, r_state, r_step, r_set_initial_state)
}

dust_cpu_seirhdage_state <- function(ptr, r_index) {
  .Call(`_seir_dust_cpu_seirhdage_state`, ptr, r_index)
}

dust_cpu_seirhdage_step <- function(ptr) {
  .Call(`_seir_dust_cpu_seirhdage_step`, ptr)
}

dust_cpu_seirhdage_reorder <- function(ptr, r_index) {
  invisible(.Call(`_seir_dust_cpu_seirhdage_reorder`, ptr, r_index))
}

dust_cpu_seirhdage_resample <- function(ptr, r_weights) {
  .Call(`_seir_dust_cpu_seirhdage_resample`, ptr, r_weights)
}

dust_cpu_seirhdage_rng_state <- function(ptr, first_only, last_only) {
  .Call(`_seir_dust_cpu_seirhdage_rng_state`, ptr, first_only, last_only)
}

dust_cpu_seirhdage_set_rng_state <- function(ptr, rng_state) {
  .Call(`_seir_dust_cpu_seirhdage_set_rng_state`, ptr, rng_state)
}

dust_cpu_seirhdage_set_data <- function(ptr, data) {
  .Call(`_seir_dust_cpu_seirhdage_set_data`, ptr, data)
}

dust_cpu_seirhdage_compare_data <- function(ptr) {
  .Call(`_seir_dust_cpu_seirhdage_compare_data`, ptr)
}

dust_cpu_seirhdage_filter <- function(ptr, step_end, save_trajectories, step_snapshot, min_log_likelihood) {
  .Call(`_seir_dust_cpu_seirhdage_filter`, ptr, step_end, save_trajectories, step_snapshot, min_log_likelihood)
}

dust_cpu_seirhdage_set_n_threads <- function(ptr, n_threads) {
  invisible(.Call(`_seir_dust_cpu_seirhdage_set_n_threads`, ptr, n_threads))
}

dust_cpu_seirhdage_n_state <- function(ptr) {
  .Call(`_seir_dust_cpu_seirhdage_n_state`, ptr)
}

dust_seirhdagevax_capabilities <- function() {
  .Call(`_seir_dust_seirhdagevax_capabilities`)
}

dust_seirhdagevax_gpu_info <- function() {
  .Call(`_seir_dust_seirhdagevax_gpu_info`)
}

dust_cpu_seirhdagevax_alloc <- function(r_pars, pars_multi, step, r_n_particles, n_threads, r_seed, deterministic, gpu_config) {
  .Call(`_seir_dust_cpu_seirhdagevax_alloc`, r_pars, pars_multi, step, r_n_particles, n_threads, r_seed, deterministic, gpu_config)
}

dust_cpu_seirhdagevax_run <- function(ptr, step_end) {
  .Call(`_seir_dust_cpu_seirhdagevax_run`, ptr, step_end)
}

dust_cpu_seirhdagevax_simulate <- function(ptr, step_end) {
  .Call(`_seir_dust_cpu_seirhdagevax_simulate`, ptr, step_end)
}

dust_cpu_seirhdagevax_set_index <- function(ptr, r_index) {
  .Call(`_seir_dust_cpu_seirhdagevax_set_index`, ptr, r_index)
}

dust_cpu_seirhdagevax_update_state <- function(ptr, r_pars, r_state, r_step, r_set_initial_state) {
  .Call(`_seir_dust_cpu_seirhdagevax_update_state`, ptr, r_pars, r_state, r_step, r_set_initial_state)
}

dust_cpu_seirhdagevax_state <- function(ptr, r_index) {
  .Call(`_seir_dust_cpu_seirhdagevax_state`, ptr, r_index)
}

dust_cpu_seirhdagevax_step <- function(ptr) {
  .Call(`_seir_dust_cpu_seirhdagevax_step`, ptr)
}

dust_cpu_seirhdagevax_reorder <- function(ptr, r_index) {
  invisible(.Call(`_seir_dust_cpu_seirhdagevax_reorder`, ptr, r_index))
}

dust_cpu_seirhdagevax_resample <- function(ptr, r_weights) {
  .Call(`_seir_dust_cpu_seirhdagevax_resample`, ptr, r_weights)
}

dust_cpu_seirhdagevax_rng_state <- function(ptr, first_only, last_only) {
  .Call(`_seir_dust_cpu_seirhdagevax_rng_state`, ptr, first_only, last_only)
}

dust_cpu_seirhdagevax_set_rng_state <- function(ptr, rng_state) {
  .Call(`_seir_dust_cpu_seirhdagevax_set_rng_state`, ptr, rng_state)
}

dust_cpu_seirhdagevax_set_data <- function(ptr, data) {
  .Call(`_seir_dust_cpu_seirhdagevax_set_data`, ptr, data)
}

dust_cpu_seirhdagevax_compare_data <- function(ptr) {
  .Call(`_seir_dust_cpu_seirhdagevax_compare_data`, ptr)
}

dust_cpu_seirhdagevax_filter <- function(ptr, step_end, save_trajectories, step_snapshot, min_log_likelihood) {
  .Call(`_seir_dust_cpu_seirhdagevax_filter`, ptr, step_end, save_trajectories, step_snapshot, min_log_likelihood)
}

dust_cpu_seirhdagevax_set_n_threads <- function(ptr, n_threads) {
  invisible(.Call(`_seir_dust_cpu_seirhdagevax_set_n_threads`, ptr, n_threads))
}

dust_cpu_seirhdagevax_n_state <- function(ptr) {
  .Call(`_seir_dust_cpu_seirhdagevax_n_state`, ptr)
}

dust_seirhdagevaxmultistrain_capabilities <- function() {
  .Call(`_seir_dust_seirhdagevaxmultistrain_capabilities`)
}

dust_seirhdagevaxmultistrain_gpu_info <- function() {
  .Call(`_seir_dust_seirhdagevaxmultistrain_gpu_info`)
}

dust_cpu_seirhdagevaxmultistrain_alloc <- function(r_pars, pars_multi, step, r_n_particles, n_threads, r_seed, deterministic, gpu_config) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrain_alloc`, r_pars, pars_multi, step, r_n_particles, n_threads, r_seed, deterministic, gpu_config)
}

dust_cpu_seirhdagevaxmultistrain_run <- function(ptr, step_end) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrain_run`, ptr, step_end)
}

dust_cpu_seirhdagevaxmultistrain_simulate <- function(ptr, step_end) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrain_simulate`, ptr, step_end)
}

dust_cpu_seirhdagevaxmultistrain_set_index <- function(ptr, r_index) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrain_set_index`, ptr, r_index)
}

dust_cpu_seirhdagevaxmultistrain_update_state <- function(ptr, r_pars, r_state, r_step, r_set_initial_state) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrain_update_state`, ptr, r_pars, r_state, r_step, r_set_initial_state)
}

dust_cpu_seirhdagevaxmultistrain_state <- function(ptr, r_index) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrain_state`, ptr, r_index)
}

dust_cpu_seirhdagevaxmultistrain_step <- function(ptr) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrain_step`, ptr)
}

dust_cpu_seirhdagevaxmultistrain_reorder <- function(ptr, r_index) {
  invisible(.Call(`_seir_dust_cpu_seirhdagevaxmultistrain_reorder`, ptr, r_index))
}

dust_cpu_seirhdagevaxmultistrain_resample <- function(ptr, r_weights) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrain_resample`, ptr, r_weights)
}

dust_cpu_seirhdagevaxmultistrain_rng_state <- function(ptr, first_only, last_only) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrain_rng_state`, ptr, first_only, last_only)
}

dust_cpu_seirhdagevaxmultistrain_set_rng_state <- function(ptr, rng_state) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrain_set_rng_state`, ptr, rng_state)
}

dust_cpu_seirhdagevaxmultistrain_set_data <- function(ptr, data) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrain_set_data`, ptr, data)
}

dust_cpu_seirhdagevaxmultistrain_compare_data <- function(ptr) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrain_compare_data`, ptr)
}

dust_cpu_seirhdagevaxmultistrain_filter <- function(ptr, step_end, save_trajectories, step_snapshot, min_log_likelihood) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrain_filter`, ptr, step_end, save_trajectories, step_snapshot, min_log_likelihood)
}

dust_cpu_seirhdagevaxmultistrain_set_n_threads <- function(ptr, n_threads) {
  invisible(.Call(`_seir_dust_cpu_seirhdagevaxmultistrain_set_n_threads`, ptr, n_threads))
}

dust_cpu_seirhdagevaxmultistrain_n_state <- function(ptr) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrain_n_state`, ptr)
}

dust_seirhdagevaxmultistrainsero_capabilities <- function() {
  .Call(`_seir_dust_seirhdagevaxmultistrainsero_capabilities`)
}

dust_seirhdagevaxmultistrainsero_gpu_info <- function() {
  .Call(`_seir_dust_seirhdagevaxmultistrainsero_gpu_info`)
}

dust_cpu_seirhdagevaxmultistrainsero_alloc <- function(r_pars, pars_multi, step, r_n_particles, n_threads, r_seed, deterministic, gpu_config) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrainsero_alloc`, r_pars, pars_multi, step, r_n_particles, n_threads, r_seed, deterministic, gpu_config)
}

dust_cpu_seirhdagevaxmultistrainsero_run <- function(ptr, step_end) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrainsero_run`, ptr, step_end)
}

dust_cpu_seirhdagevaxmultistrainsero_simulate <- function(ptr, step_end) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrainsero_simulate`, ptr, step_end)
}

dust_cpu_seirhdagevaxmultistrainsero_set_index <- function(ptr, r_index) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrainsero_set_index`, ptr, r_index)
}

dust_cpu_seirhdagevaxmultistrainsero_update_state <- function(ptr, r_pars, r_state, r_step, r_set_initial_state) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrainsero_update_state`, ptr, r_pars, r_state, r_step, r_set_initial_state)
}

dust_cpu_seirhdagevaxmultistrainsero_state <- function(ptr, r_index) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrainsero_state`, ptr, r_index)
}

dust_cpu_seirhdagevaxmultistrainsero_step <- function(ptr) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrainsero_step`, ptr)
}

dust_cpu_seirhdagevaxmultistrainsero_reorder <- function(ptr, r_index) {
  invisible(.Call(`_seir_dust_cpu_seirhdagevaxmultistrainsero_reorder`, ptr, r_index))
}

dust_cpu_seirhdagevaxmultistrainsero_resample <- function(ptr, r_weights) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrainsero_resample`, ptr, r_weights)
}

dust_cpu_seirhdagevaxmultistrainsero_rng_state <- function(ptr, first_only, last_only) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrainsero_rng_state`, ptr, first_only, last_only)
}

dust_cpu_seirhdagevaxmultistrainsero_set_rng_state <- function(ptr, rng_state) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrainsero_set_rng_state`, ptr, rng_state)
}

dust_cpu_seirhdagevaxmultistrainsero_set_data <- function(ptr, data) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrainsero_set_data`, ptr, data)
}

dust_cpu_seirhdagevaxmultistrainsero_compare_data <- function(ptr) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrainsero_compare_data`, ptr)
}

dust_cpu_seirhdagevaxmultistrainsero_filter <- function(ptr, step_end, save_trajectories, step_snapshot, min_log_likelihood) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrainsero_filter`, ptr, step_end, save_trajectories, step_snapshot, min_log_likelihood)
}

dust_cpu_seirhdagevaxmultistrainsero_set_n_threads <- function(ptr, n_threads) {
  invisible(.Call(`_seir_dust_cpu_seirhdagevaxmultistrainsero_set_n_threads`, ptr, n_threads))
}

dust_cpu_seirhdagevaxmultistrainsero_n_state <- function(ptr) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrainsero_n_state`, ptr)
}

dust_seirhdagevaxmultistrainserotimedepbeta_capabilities <- function() {
  .Call(`_seir_dust_seirhdagevaxmultistrainserotimedepbeta_capabilities`)
}

dust_seirhdagevaxmultistrainserotimedepbeta_gpu_info <- function() {
  .Call(`_seir_dust_seirhdagevaxmultistrainserotimedepbeta_gpu_info`)
}

dust_cpu_seirhdagevaxmultistrainserotimedepbeta_alloc <- function(r_pars, pars_multi, step, r_n_particles, n_threads, r_seed, deterministic, gpu_config) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrainserotimedepbeta_alloc`, r_pars, pars_multi, step, r_n_particles, n_threads, r_seed, deterministic, gpu_config)
}

dust_cpu_seirhdagevaxmultistrainserotimedepbeta_run <- function(ptr, step_end) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrainserotimedepbeta_run`, ptr, step_end)
}

dust_cpu_seirhdagevaxmultistrainserotimedepbeta_simulate <- function(ptr, step_end) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrainserotimedepbeta_simulate`, ptr, step_end)
}

dust_cpu_seirhdagevaxmultistrainserotimedepbeta_set_index <- function(ptr, r_index) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrainserotimedepbeta_set_index`, ptr, r_index)
}

dust_cpu_seirhdagevaxmultistrainserotimedepbeta_update_state <- function(ptr, r_pars, r_state, r_step, r_set_initial_state) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrainserotimedepbeta_update_state`, ptr, r_pars, r_state, r_step, r_set_initial_state)
}

dust_cpu_seirhdagevaxmultistrainserotimedepbeta_state <- function(ptr, r_index) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrainserotimedepbeta_state`, ptr, r_index)
}

dust_cpu_seirhdagevaxmultistrainserotimedepbeta_step <- function(ptr) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrainserotimedepbeta_step`, ptr)
}

dust_cpu_seirhdagevaxmultistrainserotimedepbeta_reorder <- function(ptr, r_index) {
  invisible(.Call(`_seir_dust_cpu_seirhdagevaxmultistrainserotimedepbeta_reorder`, ptr, r_index))
}

dust_cpu_seirhdagevaxmultistrainserotimedepbeta_resample <- function(ptr, r_weights) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrainserotimedepbeta_resample`, ptr, r_weights)
}

dust_cpu_seirhdagevaxmultistrainserotimedepbeta_rng_state <- function(ptr, first_only, last_only) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrainserotimedepbeta_rng_state`, ptr, first_only, last_only)
}

dust_cpu_seirhdagevaxmultistrainserotimedepbeta_set_rng_state <- function(ptr, rng_state) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrainserotimedepbeta_set_rng_state`, ptr, rng_state)
}

dust_cpu_seirhdagevaxmultistrainserotimedepbeta_set_data <- function(ptr, data) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrainserotimedepbeta_set_data`, ptr, data)
}

dust_cpu_seirhdagevaxmultistrainserotimedepbeta_compare_data <- function(ptr) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrainserotimedepbeta_compare_data`, ptr)
}

dust_cpu_seirhdagevaxmultistrainserotimedepbeta_filter <- function(ptr, step_end, save_trajectories, step_snapshot, min_log_likelihood) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrainserotimedepbeta_filter`, ptr, step_end, save_trajectories, step_snapshot, min_log_likelihood)
}

dust_cpu_seirhdagevaxmultistrainserotimedepbeta_set_n_threads <- function(ptr, n_threads) {
  invisible(.Call(`_seir_dust_cpu_seirhdagevaxmultistrainserotimedepbeta_set_n_threads`, ptr, n_threads))
}

dust_cpu_seirhdagevaxmultistrainserotimedepbeta_n_state <- function(ptr) {
  .Call(`_seir_dust_cpu_seirhdagevaxmultistrainserotimedepbeta_n_state`, ptr)
}

dust_sir_capabilities <- function() {
  .Call(`_seir_dust_sir_capabilities`)
}

dust_sir_gpu_info <- function() {
  .Call(`_seir_dust_sir_gpu_info`)
}

dust_cpu_sir_alloc <- function(r_pars, pars_multi, step, r_n_particles, n_threads, r_seed, deterministic, gpu_config) {
  .Call(`_seir_dust_cpu_sir_alloc`, r_pars, pars_multi, step, r_n_particles, n_threads, r_seed, deterministic, gpu_config)
}

dust_cpu_sir_run <- function(ptr, step_end) {
  .Call(`_seir_dust_cpu_sir_run`, ptr, step_end)
}

dust_cpu_sir_simulate <- function(ptr, step_end) {
  .Call(`_seir_dust_cpu_sir_simulate`, ptr, step_end)
}

dust_cpu_sir_set_index <- function(ptr, r_index) {
  .Call(`_seir_dust_cpu_sir_set_index`, ptr, r_index)
}

dust_cpu_sir_update_state <- function(ptr, r_pars, r_state, r_step, r_set_initial_state) {
  .Call(`_seir_dust_cpu_sir_update_state`, ptr, r_pars, r_state, r_step, r_set_initial_state)
}

dust_cpu_sir_state <- function(ptr, r_index) {
  .Call(`_seir_dust_cpu_sir_state`, ptr, r_index)
}

dust_cpu_sir_step <- function(ptr) {
  .Call(`_seir_dust_cpu_sir_step`, ptr)
}

dust_cpu_sir_reorder <- function(ptr, r_index) {
  invisible(.Call(`_seir_dust_cpu_sir_reorder`, ptr, r_index))
}

dust_cpu_sir_resample <- function(ptr, r_weights) {
  .Call(`_seir_dust_cpu_sir_resample`, ptr, r_weights)
}

dust_cpu_sir_rng_state <- function(ptr, first_only, last_only) {
  .Call(`_seir_dust_cpu_sir_rng_state`, ptr, first_only, last_only)
}

dust_cpu_sir_set_rng_state <- function(ptr, rng_state) {
  .Call(`_seir_dust_cpu_sir_set_rng_state`, ptr, rng_state)
}

dust_cpu_sir_set_data <- function(ptr, data) {
  .Call(`_seir_dust_cpu_sir_set_data`, ptr, data)
}

dust_cpu_sir_compare_data <- function(ptr) {
  .Call(`_seir_dust_cpu_sir_compare_data`, ptr)
}

dust_cpu_sir_filter <- function(ptr, step_end, save_trajectories, step_snapshot, min_log_likelihood) {
  .Call(`_seir_dust_cpu_sir_filter`, ptr, step_end, save_trajectories, step_snapshot, min_log_likelihood)
}

dust_cpu_sir_set_n_threads <- function(ptr, n_threads) {
  invisible(.Call(`_seir_dust_cpu_sir_set_n_threads`, ptr, n_threads))
}

dust_cpu_sir_n_state <- function(ptr) {
  .Call(`_seir_dust_cpu_sir_n_state`, ptr)
}
