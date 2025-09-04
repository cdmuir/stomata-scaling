# Derive V* for H1
source("r/header.R")

read_rds("objects/parameters.rds") |>
  list2env(envir = .GlobalEnv)

# this shows that there is plane of equilibria r_xy = -1 and cov(x,y) / sigma_x^2 = - 1 / beta
get_V_star(
  x = c(100, 0, 100),
  G = matrix(c(0.002, -0.001, -0.001, 0.002), 2, 2),
  V_gopt = macro_pars$Sigma[2,2] / (2 * sqrt(macro_pars$Alpha[2,2])),
  beta = constraints$beta,
  omega = micro_pars$omega,
  model = "h1"
)

