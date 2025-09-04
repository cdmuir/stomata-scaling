# Derive V* for H2
source("r/header.R")

read_rds("objects/parameters.rds") |>
  list2env(envir = .GlobalEnv)

V_gopt = macro_pars$Sigma[2, 2] / (2 * sqrt(macro_pars$Alpha[2, 2]))

G = get_msb(
  micro_pars$phi_f,
  micro_pars$omega,
  micro_pars$M,
  constraints$beta,
  x = 1e2 * micro_pars$M,
  model = "h2"
)

beta = 0.5
omega = micro_pars$omega
phi_f = micro_pars$phi_f

V0 = matrix(c(4 * V_gopt, rep(-4 * V_gopt + 1e-5, 2), 4 * V_gopt), 2, 2)

for (i in 1:10000) {
  V0 = matrix(rep(c(
    get_Vd_h2(V0, V_gopt, G, beta, omega, phi_f),
    get_Vda_h2(V0, V_gopt, G, beta, omega, phi_f),
    get_Va_h2(V0, V_gopt, G, beta, omega, phi_f)
  ), c(1, 2, 1)), 2, 2)
}

V0

get_V_star(
  x = c(1.6, -1.6, 1.6),
  G = G,
  V_gopt = macro_pars$Sigma[2, 2] / (2 * sqrt(macro_pars$Alpha[2, 2])),
  beta = constraints$beta,
  omega = micro_pars$omega,
  phi_f = micro_pars$phi_f,
  model = "h2"
)
