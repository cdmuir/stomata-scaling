# For each fitted simulation, pull parameter estimates for parametric bootstrap CIs
source("r/header.R")

read_rds("objects/parameters.rds") |>
  list2env(envir = .GlobalEnv)

df_ouch = list.files("objects/ouch") |>
  map_dfr(\(.x) {
    fit_ouch = read_rds(paste0("objects/ouch/", .x))
    
    pars = fit_ouch$FinalFound$ParamsInModel
    V = fit_ouch$FinalFound$ParamSummary$stationary.cov.matrix
    
    tibble(
      sim = .x,
      A = list(pars$A),
      mPsi = list(pars$mPsi),
      vY0 = list(pars$vY0),
      Syy = list(pars$Syy),
      V = list(V)
    )
    
  })

write_rds(df_ouch, "objects/df_ouch.rds")

# Values for Table 3
tibble(
  dS_angio = map_dbl(df_ouch$mPsi, extract2, 1, 1),
  dS_gymno = map_dbl(df_ouch$mPsi, extract2, 1, 2),
  dS_pterido = map_dbl(df_ouch$mPsi, extract2, 1, 3),
  aS_angio = map_dbl(df_ouch$mPsi, extract2, 2, 1),
  aS_gymno = map_dbl(df_ouch$mPsi, extract2, 2, 2),
  aS_pterido = map_dbl(df_ouch$mPsi, extract2, 2, 3)
) |>
  pivot_longer(cols = everything(), names_to = c("trait", "group"), names_sep = "_", values_to = "value") |>
  group_by(trait, group) |>
  point_interval(value) |>
  mutate(across(value:.upper, exp))

tibble(
  dS = map_dbl(df_ouch$V, extract2, 1, 1),
  aS = map_dbl(df_ouch$V, extract2, 2, 2),
  dSaS = map_dbl(df_ouch$V, extract2, 1, 2),
  sma = -sign(dSaS) * sqrt(aS / dS),
  ols1 = -dS / dSaS,
  ols2 = -dSaS / aS,
  gopt = constraints$beta ^ 2 * aS + dS + 2 * constraints$beta * dSaS,
  aSgopt = constraints$beta * aS + dSaS
) |>
  pivot_longer(cols = everything()) |>
  group_by(name) |>
  point_interval(value) 

