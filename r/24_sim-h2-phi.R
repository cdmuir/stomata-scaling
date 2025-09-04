# Recursion equations to determine stationary variance of V* under with H2 variable phi_f
source("r/header.R")

read_rds("objects/parameters.rds") |>
  list2env(envir = .GlobalEnv)

# Start from a single common ancestor and simulate trait evolution on a star
# phylogeny with n_spp species.

dir1 = "sims/h2_phi"
if (!dir.exists(dir1)) {
  dir.create(dir1)
}

n_spp = n_spp
set.seed(96780258)
seeds = sample(1e9, n_spp)

# Set 1. phi_f = 0.5 * micro_pars$phi_f
# Set 2. phi_f = 2.0 * micro_pars$phi_f

mip = replicate(2, micro_pars, simplify = FALSE)
mip[[1]]$phi_f = 0.5 * micro_pars$phi_f
mip[[2]]$phi_f = 2 * micro_pars$phi_f

plan(multisession, workers = n_workers)

for (s in seq_along(mip)) {
  future_walk(
    seq_len(n_spp),
    \(.x, dir1) {
      n = str_pad(.x, ceiling(log10(n_spp) + 1), pad = 0)
      
      simulate_lineage(
        n_gen = 5e6,
        constraints = constraints,
        micro_pars = mip[[s]],
        macro_pars = macro_pars,
        x = 1e1 * mip[[s]]$M,
        model = "h2",
        individual_based = FALSE,
        drift = FALSE,
        thin = 5e3,
        init_method = "Theta",
        custom_init = NULL,
        n_ind_sim = 1,
        seed = seeds[.x] + s
      ) |>
        mutate(n = n, phi_f = mip[[s]]$phi_f) |>
        write_rds(glue(
          "{dir1}/sim_{set}_{n}.rds",
          set = paste0("M", s),
          n = n
        ))
    },
    dir1 = dir1,
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
  )
  
}

files_to_zip = list.files(dir1, full.names = TRUE)
zip(zipfile = str_c(dir1, ".zip"), files = files_to_zip)
