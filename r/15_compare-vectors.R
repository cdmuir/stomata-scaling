source("r/header.R")

read_rds("objects/parameters.rds") |>
  list2env(envir = .GlobalEnv)

rM_names = c("Negative correlation",
             "No correlation",
             "Positive correlation")
MdS_names = str_c(c("Low", "Medium", "High"), " $M_{d_\\mathrm{S}}$")
MaS_names = str_c(c("Low", "Medium", "High"), " $M_{a_\\mathrm{S}}$")
phif_names = str_c(c("Weak", "Moderate", "Strong"), " selection on $f_\\textrm{S}$")

plan(multisession, workers = 10)
df1 = crossing(
  nesting(
    M_dS = micro_pars$M[1, 1] * 10^seq(-1, 1, by = 1),
    MdS_names = factor(MdS_names, levels = MdS_names)
  ),
  nesting(
    M_aS = micro_pars$M[2, 2] * 10^seq(-1, 1, by = 1),
    MaS_names = factor(MaS_names, levels = MaS_names)
  ),
  nesting(
    r_M = seq(-0.5, 0.5, 0.5),
    rM_names = factor(rM_names, levels = rM_names)
  ),
  omega = micro_pars$omega,
  phi_a = 10^seq(-1, 1, 1),
  nesting(
    phi_f = 10^seq(-1, 1, 1),
    phif_names = factor(phif_names, levels = phif_names)
  ),
  beta = constraints$beta,
) |>
  mutate(M_dSaS = r_M * sqrt(M_dS * M_aS)) |>
  rowwise() |>
  mutate(M = list(matrix(c(
    M_dS, M_dSaS, M_dSaS, M_aS
  ), 2, 2, byrow = TRUE)), ) %>%
  split(~ seq_len(nrow(.))) |>
  future_map_dfr(\(.x) {
    G_h1 = matrix(c(1, rep(-0.5, 2), 0.25), 2, 2)
    G_h2 = get_Gstar(.x$phi_f,
                   .x$omega,
                   .x$M[[1]],
                   constraints$beta,
                   model = "h2")
    G_h3 = get_Gstar(.x$phi_a,
                   .x$omega,
                   .x$M[[1]],
                   constraints$beta,
                   model = "h3")
    
    # Starting trait values
    gmax_h1 = macro_pars$Theta[2, 1]
    da_bar_h1 = c(
      gmax_h1 - constraints$z0_gmax - constraints$beta * macro_pars$Theta[1, 1],
      macro_pars$Theta[1, 1]
    )
    
    gmax_h2 = macro_pars$Theta[2, 1]
    da_bar_h2 = c(
      (
        constraints$beta * constraints$fS_min - constraints$beta * constraints$z0_fS - gmax_h2 + constraints$z0_gmax
      ) / (constraints$beta - 1),
      (
        -constraints$fS_min + gmax_h2 + constraints$z0_fS - constraints$z0_gmax
      ) / (constraints$beta - 1)
    )
    
    gmax_h3 = macro_pars$Theta[2, 1]
    g_opt = gmax_h3 + 1 / (2 * .x$phi_f)
    aS = macro_pars$Theta[1, 1]
    aS_opt = (aS * 2 * .x$phi_f + constraints$beta * .x$phi_a + .x$phi_a) / (2 * .x$phi_f)
    
    da_bar_h3 = c(
      -aS_opt * constraints$beta - constraints$beta ** 2 * .x$phi_a / (2 * .x$phi_f) + constraints$beta * .x$phi_a / (2 * .x$phi_f) + g_opt - constraints$z0_gmax - 1 / (2 * .x$phi_f),
      aS
    )
    
    d_gopt = sdOU(1, sqrt(macro_pars$Alpha[2, 2]), sqrt(macro_pars$Sigma[2, 2]))
    
    del_logW1_h1 = get_del_logW(
      gmax_h1 + d_gopt,
      da_bar_h1,
      .x$omega,
      constraints$z0_gmax,
      constraints$beta,
      model = "h1"
    )
    d_h1 = G_h1 %*% del_logW1_h1
    
    del_logW1_h2 = get_del_logW(
      gmax_h2 + d_gopt,
      constraints$fS_min,
      da_bar_h2,
      .x$phi_f,
      .x$omega,
      constraints$z0_fS,
      constraints$z0_gmax,
      constraints$beta,
      model = "h2"
    )
    
    d_h2 = G_h2 %*% del_logW1_h2
    
    del_logW1_h3 = get_del_logW(
      c(aS_opt, gmax_h3 + d_gopt),
      da_bar_h3,
      .x$phi_a,
      .x$phi_f,
      .x$omega,
      constraints$z0_gmax,
      constraints$beta,
      model = "h3"
    )
    
    d_h3 = G_h3 %*% del_logW1_h3
    
    tibble(
      M_dS = .x$M_dS,
      MdS_names = .x$MdS_names,
      M_aS = .x$M_aS,
      MaS_names = .x$MaS_names,
      r_M = .x$r_M,
      rM_names = .x$rM_names,
      omega = .x$omega,
      phi_a = .x$phi_a,
      phi_f = .x$phi_f,
      phif_names = .x$phif_names,
      h = c("h1", "h2", "h3"),
      hypothesis = factor(hypotheses, levels = hypotheses),
      d = c(d_h1[2] / d_h1[1], d_h2[2] / d_h2[1], d_h3[2] / d_h3[1])
    )
    
  }, .progress = TRUE) |>
  filter(!is.na(d)) |>
  mutate(angle = atan(d))

# Effect of mutational correlation
fig_vec1 = df1 |>
  filter(
    phi_a == micro_pars$phi_a,
    phi_f == micro_pars$phi_f,
    M_dS == micro_pars$M[1, 1],
    M_aS == micro_pars$M[2, 2]
  ) |>
  ggplot(aes(
    x = 0,
    y = 0,
    angle = angle,
    radius = 1,
    color = hypothesis
  )) +
  facet_wrap(~ rM_names) +
  geom_spoke(arrow = arrow(length = unit(0.25, "cm")), linewidth = 1.1) +
  xlab("Relative response of stomatal density") +
  ylab("Relative response of stomatal size") +
  scale_x_continuous(breaks = seq(0, 1, 0.5)) +
  ylim(-1, 1) +
  coord_equal() +
  theme(
    legend.position = "bottom",
    legend.direction = "vertical",
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Effect of mutational asymmetry
fig_vec2 = df1 |>
  filter(r_M == 0, phi_a == micro_pars$phi_a, phi_f == micro_pars$phi_f, ) |>
  ggplot(aes(
    x = 0,
    y = 0,
    angle = angle,
    radius = 1,
    color = hypothesis
  )) +
  facet_grid(MdS_names ~ MaS_names) +
  geom_spoke(arrow = arrow(length = unit(0.25, "cm")), linewidth = 1.1) +
  xlab("Relative response of stomatal density") +
  ylab("Relative response of stomatal size") +
  scale_x_continuous(breaks = seq(0, 1, 0.5)) +
  ylim(-1, 1) +
  coord_equal() +
  theme(
    legend.position = "bottom",
    legend.direction = "vertical",
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Effect of phi_f
fig_vec3 = df1 |>
  filter(r_M == 0,
         phi_a == micro_pars$phi_a,
         M_dS == micro_pars$M[1, 1],
         M_aS == micro_pars$M[2, 2]) |>
  ggplot(aes(
    x = 0,
    y = 0,
    angle = angle,
    radius = 1,
    color = hypothesis
  )) +
  facet_wrap(~ phif_names) +
  geom_spoke(arrow = arrow(length = unit(0.25, "cm")), linewidth = 1.1) +
  xlab("Relative response of stomatal density") +
  ylab("Relative response of stomatal size") +
  scale_x_continuous(breaks = seq(0, 1, 0.5)) +
  ylim(-1, 1) +
  coord_equal() +
  theme(
    legend.position = "bottom",
    legend.direction = "vertical",
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Render figures
tikz(file = "ms/figures/fig_vec1.tex",
     width = 6.5,
     height = 6)
print(fig_vec1)
dev.off()

tikz(file = "ms/figures/fig_vec2.tex",
     width = 6.5,
     height = 8)
print(fig_vec2)
dev.off()

tikz(file = "ms/figures/fig_vec3.tex",
     width = 6.5,
     height = 6)
print(fig_vec3)
dev.off()
