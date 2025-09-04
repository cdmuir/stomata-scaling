# Recursion equations to determine stationary variance of V* under with H2 variable rho_M
# n.b. - individual simulations are not included with source code
source("r/header.R")

read_rds("objects/parameters.rds") |>
  list2env(envir = .GlobalEnv)

# Start from a single common ancestor and simulate trait evolution on a star
# phylogeny with n_spp species.

dir1 = "sims/h2_rM"
if (!dir.exists(dir1)) {
  dir.create(dir1)
}

n_spp = n_spp #1e4
set.seed(628682545)
seeds = sample(1e9, n_spp)

# Set 1. rho_M = -0.5
# Set 2. rho_M = 0.5

mip = replicate(2, micro_pars, simplify = FALSE)
mip[[1]]$M[1, 2] = mip[[1]]$M[2, 1] = -0.5 * sqrt(prod(diag(micro_pars$M)))
mip[[2]]$M[1, 2] = mip[[2]]$M[2, 1] = 0.5 * sqrt(prod(diag(micro_pars$M)))

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
        x = 1e2 * mip[[s]]$M,
        model = "h2",
        individual_based = FALSE,
        drift = FALSE,
        thin = 5e3,
        init_method = "Theta",
        custom_init = NULL,
        n_ind_sim = 1,
        seed = seeds[.x] + s
      ) |>
        mutate(n = n, M = list(mip[[s]]$M)) |>
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
