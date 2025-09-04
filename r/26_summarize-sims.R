source("r/header.R")

set.seed(423028207)

read_rds("objects/parameters.rds") |>
  list2env(envir = .GlobalEnv)

# Load simulation output
sets = list.files("sims", pattern = "^h[23]_")

walk(sets, \(.s) {
  assign(.s,
         list.files(paste0("sims/", .s), full.names = TRUE) |>
           map_dfr(read_rds),
         envir = .GlobalEnv)
})

# Baseline scenarios ----
df_baseline = bind_rows(
  h2_baseline |>
    filter(gen == max(gen)) |>
    select(n, matches("^dag_bar_[1-3]$")) |>
    mutate(hypothesis = "h2", phi_a = NA),
  h3_baseline |>
    filter(gen == max(gen)) |>
    select(n, matches("^dag_bar_[1-3]$")) |>
    mutate(hypothesis = "h3", phi_a = micro_pars$phi_a)
) |>
  mutate(
    M_dS = micro_pars$M[1, 1],
    M_aS = micro_pars$M[2, 2],
    Mratio = M_dS / M_aS,
    rho_M = 0,
    omega = micro_pars$omega,
    phi_f = micro_pars$phi_f,
  )

# 1. M-ratio ----

bind_rows(
  h2_Mratio |>
    filter(gen == max(gen)) |>
    select(n, matches("^dag_bar_[1-3]$"), M) |>
    mutate(hypothesis = "h2"),
  h3_Mratio |>
    filter(gen == max(gen)) |>
    select(n, matches("^dag_bar_[1-3]$"), M) |>
    mutate(hypothesis = "h3")
) |>
  rowwise() |>
  mutate(M_dS = M[1, 1],
         M_aS = M[2, 2],
         Mratio = M_dS / M_aS) |>
  bind_rows(select(
    df_baseline,
    n,
    matches("^dag_bar_[1-3]$"),
    M_dS,
    M_aS,
    Mratio,
    hypothesis
  )) |>
  ungroup() |>
  rename(
    dS = dag_bar_1,
    aS = dag_bar_2,
    gmax = dag_bar_3
  ) |>
  summarize(
    boot_var(dS),
    boot_var(aS),
    boot_var(gmax),
    boot_cov(dS, aS),
    .by = c("hypothesis", "M_dS", "M_aS", "Mratio")
  ) |>
  pivot_longer(
    cols = matches("^V_"),
    names_pattern = "^(V_[A-Za-z]+)_([a-z]+)$",
    names_to = c("component", "quantity"),
    values_to = "value"
  ) |>
  write_rds("objects/df_Mratio.rds")

# 2. rho_M ----

bind_rows(
  h2_rM |>
    filter(gen == max(gen)) |>
    select(n, matches("^dag_bar_[1-3]$"), M) |>
    mutate(hypothesis = "h2"),
  h3_rM |>
    filter(gen == max(gen)) |>
    select(n, matches("^dag_bar_[1-3]$"), M) |>
    mutate(hypothesis = "h3"),
  df_baseline |> 
    select(n, matches("^dag_bar_[1-3]$"), hypothesis) |>
    mutate(M = list(micro_pars$M))
) |>
  rowwise() |>
  mutate(rho_M = M[1, 2] / sqrt(prod(M[1, 1], M[2, 2]))) |>
  ungroup() |>
  rename(
    dS = dag_bar_1,
    aS = dag_bar_2,
    gmax = dag_bar_3
  ) |>
  summarize(
    boot_var(dS),
    boot_var(aS),
    boot_var(gmax),
    boot_cov(dS, aS),
    .by = c("hypothesis", "rho_M")
  ) |>
  pivot_longer(
    cols = matches("^V_"),
    names_pattern = "^(V_[A-Za-z]+)_([a-z]+)$",
    names_to = c("component", "quantity"),
    values_to = "value"
  ) |>
  write_rds("objects/df_rM.rds")

# 3. omega ----
  
bind_rows(
  h2_omega |>
    filter(gen == max(gen)) |>
    select(n, matches("^dag_bar_[1-3]$"), omega) |>
    mutate(hypothesis = "h2"),
  h3_omega |>
    filter(gen == max(gen)) |>
    select(n, matches("^dag_bar_[1-3]$"), omega) |>
    mutate(hypothesis = "h3"),
  df_baseline |> select(n, matches("^dag_bar_[1-3]$"), hypothesis, omega)
) |>
  rename(
    dS = dag_bar_1,
    aS = dag_bar_2,
    gmax = dag_bar_3
  ) |>
  summarize(
    boot_var(dS),
    boot_var(aS),
    boot_var(gmax),
    boot_cov(dS, aS),
    .by = c("hypothesis", "omega")
  ) |>
  pivot_longer(
    cols = matches("^V_"),
    names_pattern = "^(V_[A-Za-z]+)_([a-z]+)$",
    names_to = c("component", "quantity"),
    values_to = "value"
  ) |>
  write_rds("objects/df_omega.rds")

# 4. phi ----

bind_rows(
  h2_phi |>
    filter(gen == max(gen)) |>
    select(n, matches("^dag_bar_[1-3]$"), phi_f) |>
    mutate(hypothesis = "h2", phi_a = NA),
  h3_phi |>
    filter(gen == max(gen)) |>
    select(n, matches("^dag_bar_[1-3]$"), phi_f, phi_a) |>
    mutate(hypothesis = "h3"),
  df_baseline |> select(n, matches("^dag_bar_[1-3]$"), hypothesis, phi_f, phi_a)
) |>
  rename(
    dS = dag_bar_1,
    aS = dag_bar_2,
    gmax = dag_bar_3
  ) |>
  summarize(
    boot_var(dS),
    boot_var(aS),
    boot_var(gmax),
    boot_cov(dS, aS),
    .by = c("hypothesis", "phi_f", "phi_a")
  ) |>
  pivot_longer(
    cols = matches("^V_"),
    names_pattern = "^(V_[A-Za-z]+)_([a-z]+)$",
    names_to = c("component", "quantity"),
    values_to = "value"
  ) |>
  write_rds("objects/df_phi.rds")
