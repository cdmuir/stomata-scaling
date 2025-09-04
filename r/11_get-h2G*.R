# Plots of affect on (co)variance of G* in H2.
source("r/header.R")

read_rds("objects/parameters.rds") |>
  list2env(envir = .GlobalEnv)

G_h2_pars = list(
  M_dS = micro_pars$M[1,1] * 10 ^ seq(-1, 1, by = 0.1),
  M_aS = micro_pars$M[2,2] * 10 ^ seq(-1, 1, by = 0.1),
  r_M = seq(-0.7, 0.7, 0.1),
  omega = micro_pars$omega * c(0.5, 1, 2),
  phi_f = 2 ^ seq(-1, 1, 0.25),
  beta = constraints$beta
)

df1 = crossing(
  M_dS  = G_h2_pars$M_dS,
  M_aS  = G_h2_pars$M_aS,
  r_M   = G_h2_pars$r_M,
  omega = G_h2_pars$omega,
  phi_f = G_h2_pars$phi_f,
  beta  = G_h2_pars$beta
) |>
  mutate(
    M_dSaS = r_M * sqrt(M_dS * M_aS)) |>
  rowwise() |>
  mutate(
    M = list(matrix(c(M_dS, M_dSaS, M_dSaS, M_aS), 2, 2, byrow = TRUE)),
    G_star = list(get_msb(phi_f, omega, M, beta, x = 1e2 * M, model = "h2"))
  ) |>
  filter(!is.na(G_star[1])) |>
  mutate(
    G_dS = G_star[1,1],
    G_aS = G_star[2,2],
    G_dSaS = G_star[1,2]
  )

# Effect of mutational variance and selection on trait variance
df2 = bind_rows(
  df1 |>
    filter(round(r_M, 1) == 0, phi_f == 2, M_aS == micro_pars$M[2, 2]) |>
    select(M = M_dS, G = G_dS, omega) |>
    mutate(trait = "density"),
  df1 |>
    filter(round(r_M, 1) == 0, phi_f == 2, M_dS == micro_pars$M[1, 1]) |>
    select(M = M_aS, G = G_aS, omega) |>
    mutate(trait = "size")
) |>
  mutate(`Strength of\nselection on $g_\\mathrm{s,max}$` = factor(case_when(
    omega == micro_pars$omega / 2 ~ "strong",
    omega == micro_pars$omega ~ "medium",
    omega == micro_pars$omega * 2 ~ "weak"
  ), levels = c("weak", "medium", "strong")))

figH2G1 = ggplot(df2,
                 aes(M, G, linetype = `Strength of\nselection on $g_\\mathrm{s,max}$`)) +
  facet_wrap(~ trait) +
  geom_line(linewidth = 1.1) +
  labs(x = "Mutational variance", y = "Additive genetic variance") +
  scale_x_log10(breaks = c(5e-7, 2e-6, 8e-6, 32e-6),
                labels = c("$0.5 \\times 10^{-6}$", "$2.0 \\times 10^{-6}$", "$8.0 \\times 10^{-6}$", "$3.2 \\times 10^{-5}$")) +
  scale_y_log10(limits = c(4e-4, 4e-3), breaks = c(5e-4, 1e-3, 2e-3, 4e-3),
                labels = c("$0.5 \\times 10^{-3}$", "$1.0 \\times 10^{-3}$", "$2.0 \\times 10^{-3}$", "$4.0 \\times 10^{-3}$")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Effect of mutational correlation and phi_f on trait correlation
figH2G2 = df1 |>
  filter(
    M_dS == micro_pars$M[1, 1],
    M_aS == micro_pars$M[2, 2],
    phi_f %in% c(0.5, 1, 2)
  )  |>
  mutate(
    r_G = G_dSaS / sqrt(G_dS * G_aS),
    `Strength of\nselection on\nstomatal size` = factor(case_when(
    phi_f == 0.5 ~ "strong",
    phi_f == 1 ~ "medium",
    phi_f == 2 ~ "weak"
  ), levels = c("strong", "medium", "weak"))) |>
  ggplot(aes(r_M, r_G, linetype = `Strength of\nselection on\nstomatal size`)) +
  geom_line(linewidth = 1.1) +
  xlab("Mutational correlation") +
  ylab("Additive genetic correlation") +
  xlim(-0.75, 0.75) +
  ylim(-1.0, -0.3) +
  coord_equal()

# Effect of phi_a on relative variance of density and size
M_as_levels = micro_pars$M[2,2] * 10 ^ seq(-1, 1, by = 1)
figH2G3 = df1 |>
  filter(round(r_M, 1) == 0, omega == micro_pars$omega,
         M_dS == micro_pars$M[1,1],
         M_aS %in% M_as_levels) |>
  mutate(`Mutational variance\nof stomatal size` = factor(case_when(
    M_aS == M_as_levels[1] ~ "low",
    M_aS == M_as_levels[2] ~ "medium",
    M_aS == M_as_levels[3] ~ "high"
  ), levels = c("low", "medium", "high"))) |>
  mutate(rel_var = G_dS / G_aS) |>
  ggplot(aes(phi_f, rel_var, linetype = `Mutational variance\nof stomatal size`)) +
  geom_line(linewidth = 1.1) +
  xlab("Strength of selection on $f_\\mathrm{S}$") +
  ylab("Relative genetic variance, $G^*_{d_\\mathrm{S}} / G^*_{a_\\mathrm{S}}$") +
  scale_x_log10(breaks = c(0.5, 1, 2), labels = paste0("$", c(0.5, 1, 2), "\\omega$")) +
  scale_y_log10(limits = c(1/5, 5), breaks = c(1/4, 1, 4)) +
    coord_equal()

# Write parameters and figures
tikz(file = "ms/figures/figH2G1.tex", width = 6, height = 4)
print(figH2G1)  
dev.off()

tikz(file = "ms/figures/figH2G2.tex", width = 6, height = 4)
print(figH2G2)  
dev.off()

tikz(file = "ms/figures/figH2G3.tex", width = 6, height = 4)
print(figH2G3)  
dev.off()

write_rds(G_h2_pars, "objects/G_h2_pars.rds")
