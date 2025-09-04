# Recursion equations to determine stationary variance of V* under H2 baseline scenario
# n.b. - individual simulations are not included with source code
source("r/header.R")

read_rds("objects/parameters.rds") |>
  list2env(envir = .GlobalEnv)

# Start from a single common ancestor and simulate trait evolution on a star
# phylogeny with n_spp species.

dir1 = "sims/h2_baseline1"
if (!dir.exists(dir1)) {
  dir.create(dir1)
}

n_spp = n_spp
set.seed(827897604)
seeds = sample(1e9, n_spp)

plan(multisession, workers = 10)
future_walk(
  seq_len(n_spp),
  \(.x, dir1) {
    n = str_pad(.x, ceiling(log10(n_spp) + 1), pad = 0)
    simulate_lineage(
      n_gen = 5e6,
      constraints = constraints,
      micro_pars = micro_pars,
      macro_pars = macro_pars,
      x = 1e1 * micro_pars$M,
      model = "h2",
      individual_based = FALSE,
      drift = FALSE,
      thin = 5e3,
      init_method = "Theta",
      custom_init = NULL,
      n_ind_sim = 1,
      seed = seeds[.x]
    ) |>
      mutate(n = n) |>
      write_rds(glue("{dir1}/sim_{n}.rds", n = n))
  },
  dir1 = dir1,
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)

files_to_zip = list.files(dir1, full.names = TRUE)
zip(zipfile = str_c(dir1, ".zip"), files = files_to_zip)
