# Fit synthetic data from OUCH model for parameteric bootstrap
# n.b. - individual fitted models are not included with source code

source("r/header.R")

phy = read_rds("objects/phy.rds")

fit_simmap = read_rds("objects/fit_simmap.rds")

rgs = fit_simmap$maps |>
  map_chr(\(.x) names(.x)[which.max(.x)])

plan(multisession, workers = 10)

list.files("sims/ouch") |>
  future_walk(\(.x) {
    mdat = read_rds(glue("sims/ouch/{.x}"))
    fit_ouch = ouchModel(phy, mdat, regimes = rgs, Atype = "Symmetric")
    write_rds(fit_ouch, glue("objects/ouch/{.x}"))
  }, .progress = TRUE, .options = furrr_options(seed = TRUE))
