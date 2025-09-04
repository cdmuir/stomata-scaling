# Recursion equations to determine stationary variance of V* under with H2 variable M ratio
# n.b. - individual simulations are not included with source code
source("r/header.R")

read_rds("objects/parameters.rds") |>
  list2env(envir = .GlobalEnv)

# Start from a single common ancestor and simulate trait evolution on a star
# phylogeny with n_spp species.

dir1 = "sims/h3_Mratio"
if (!dir.exists(dir1)) {
  dir.create(dir1)
}

n_spp = n_spp
set.seed(812203908)
seeds = sample(1e9, n_spp)

# Set 1. M_dS = 0.1 * M_aS
# Set 2. M_dS = 10 * M_aS
# Set 3. 0.1 * M_dS = M_aS
# Set 4. 10 * M_dS = M_aS

mip = replicate(4, micro_pars, simplify = FALSE)
mip[[1]]$M[1, 1] = 0.1 * micro_pars$M[2, 2]
mip[[2]]$M[1, 1] = 10 * micro_pars$M[2, 2]
mip[[3]]$M[2, 2] = 0.1 * micro_pars$M[1, 1]
mip[[4]]$M[2, 2] = 10 * micro_pars$M[1, 1]

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
        model = "h3",
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
