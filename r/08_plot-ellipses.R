source("r/header.R")

# Derivations for why H1 and H2 are plotted as they are. They are constrained so 
# that the variances equals the stationary variance of Vg.
# H1
# Vd = 4 Va
# r = -1
# Therefore Cov(a_s,d_s) = -2 Va
# Vg = Vd + b ^ 2 * Va + 2 * b * Cov(a_s,d_s)
# Vg = 4 * Va + b ^ 2 * Va + 2 * b * -2 Va
# Va = Vg / 2.25
# Vg = 2
# Va = Vg / 2.25
# Vd = 4 * Va
# Vad = - 2 * Va
# b = 0.5
# Vd + b ^ 2 * Va + 2 * b * Vad

# H2
# Va = Vd = 4 Vg
# r = -1
# Therefore Cov(a_s,d_s) = -4 Vg
# Vg = Vd + b ^ 2 * Va + 2 * b * Cov(a_s,d_s)
# Vg = 4 Vg + b ^ 2 * 4 Vg + 2 * b * -Vg

# Vg = 2
# Vd = 4 * Vg
# Va = 4 * Vg
# Vad = -4 * Vg
# b = 0.5
# Vd + b ^ 2 * Va + 2 * b * Vad

read_rds("objects/parameters.rds") |>
  list2env(envir = .GlobalEnv)

dat = read_rds("objects/dat.rds")

df_ouch = read_rds("objects/df_ouch.rds")

df_mu = tibble(
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
  filter(group == "angio") |>
  select(trait, value) |>
  pivot_wider(names_from = trait)

df_Sigma = tibble(
  dS = map_dbl(df_ouch$V, extract2, 1, 1),
  aS = map_dbl(df_ouch$V, extract2, 2, 2),
  dSaS = map_dbl(df_ouch$V, extract2, 1, 2),
  gopt = constraints$beta ^ 2 * aS + dS + 2 * constraints$beta * dSaS
) |>
  pivot_longer(cols = everything()) |>
  group_by(name) |>
  point_interval(value) |>
  select(name, value) |>
  pivot_wider()

df_ellipse = list(
  hypothesis = hypotheses,
  rho = c(-1, -1, -0.5),
  var_aS = c(df_Sigma$gopt / 2.25, 4 * df_Sigma$gopt, df_Sigma$aS),
  var_dS = c(4 * df_Sigma$gopt / 2.25, 4 * df_Sigma$gopt, df_Sigma$dS)
) |>
  pmap_dfr( ~ {
    ellipse(
      x = matrix(c(1, ..2, ..2, 1), ncol = 2),
      scale = c(..3, ..4),
      centre = c(df_mu$aS, df_mu$dS),
      level = 0.999
    ) |>
      as_tibble(.name_repair = NULL) |>
      rename(a_s = x, d_s = y) |>
      mutate(hypothesis = ..1,
             rho = ..2,
             var_aS = ..3,
             var_dS = ..4)
    
  }) |>
  mutate(hypothesis = factor(hypothesis, levels = hypotheses))

df_ellipse |>
  ggplot(aes(x = a_s, y = d_s)) +
  geom_point(data = dat, color = "grey", alpha = 0.5) + 
  geom_path(aes(group = hypothesis, color = hypothesis), linewidth = 1.1) +
  scale_x_continuous(
    limits = df_mu$aS + c(-10, 10) * df_Sigma$aS,
    breaks = log(c(10, 100, 1000)),
    labels = c(10, 100, 1000)
  ) + 
  scale_y_continuous(
    limits = df_mu$aS + c(-10, 10) * df_Sigma$aS,
    breaks = log(c(10, 100, 1000)),
    labels = c(10, 100, 1000)
  ) + 
  xlab(expression(paste("Stomatal size [",mu, m ^2, "], log-scale"))) +
  ylab(expression(atop(Stomatal~density, paste(group("[", pores~mm ^-2, "]"), ",",~log-scale)))) +
  ggtitle("Predicted versus observed among species trait covariance", subtitle = "Note: predicted covariance ellipses, not regression lines, are shown") +
  coord_equal() +
  theme_bw() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 12),
  )

ggsave("ms/figures/ellipses.pdf", width = 6, height = 4)  
