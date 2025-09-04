# Simulate synthetic data from OUCH model for parametric bootstrap
# n.b. - individual simulation outputs are not included with source code
source("r/header.R")

phy = read_rds("objects/phy.rds")

fit_simmap = read_rds("objects/fit_simmap.rds")
fit_ouch = read_rds("objects/fit_ouch.rds")

rgs = fit_simmap$maps |>
  map_chr(\(.x) names(.x)[which.max(.x)])

if (!dir.exists("sims/ouch"))
  dir.create("sims/ouch")

n_sim = 1e3
set.seed(420918570)

plan(multisession, workers = 10)
future_walk(seq_len(n_sim),
            \(.x) {
              sim = simulOUCHProcPhylTree(phy, fit_ouch$FinalFound$ParamsInModel, regimes = rgs)
              write_rds(sim, glue("sims/ouch/sim_{.y}.rds", .y = str_pad(.x, 4, pad = "0")))
            },
            .progress = TRUE,
            .options = furrr_options(seed = TRUE))
