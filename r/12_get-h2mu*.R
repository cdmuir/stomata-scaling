source("r/header.R")

read_rds("objects/parameters.rds") |>
  list2env(envir = .GlobalEnv)

G = get_msb(micro_pars$phi_f,
            micro_pars$omega,
            micro_pars$M,
            constraints$beta,
            x = 1e1 * micro_pars$M,
            model = "h2")

get_mu_star(
  G = G,
  omega = micro_pars$omega,
  phi_f = micro_pars$phi_f,
  z0_gmax = constraints$z0_gmax,
  z0_fS = constraints$z0_fS,
  beta = constraints$beta,
  fS_min = constraints$fS_min,
  g_opt = macro_pars$Theta[2, 1],
  model = "h2"
)

constraints$z0_fS + 8.559028 - 1.651273
constraints$z0_gmax + 8.559028 - 0.5 * 1.651273
