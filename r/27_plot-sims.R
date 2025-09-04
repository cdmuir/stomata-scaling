source("r/header.R")

read_rds("objects/parameters.rds") |>
  list2env(envir = .GlobalEnv)

read_rds("objects/hyperparameters.rds") |>
  list2env(envir = .GlobalEnv)

components = c(
  "$V^*_{a_\\mathrm{S}}$",
  "$V^*_{d_\\mathrm{S}}$",
  "$V^*_{d_\\mathrm{S},a_\\mathrm{S}}$",
  "$V^*_{g_\\mathrm{s,max}}$"
)

# Equilibria ----
Vstar_h3 = StationaryVariance(macro_pars$Alpha, macro_pars$Sigma)

df_eq = tribble(
  ~ hypothesis1,
  ~ component1,
  ~ estimate,
  "$H_2$",
  "$V^*_{d_\\mathrm{S}}$",
  macro_pars$Sigma[2, 2] / (2 * sqrt(macro_pars$Alpha[2, 2]) * constraints$beta^2),
  "$H_2$",
  "$V^*_{a_\\mathrm{S}}$",
  macro_pars$Sigma[2, 2] / (2 * sqrt(macro_pars$Alpha[2, 2]) * constraints$beta^2),
  "$H_2$",
  "$V^*_{d_\\mathrm{S},a_\\mathrm{S}}$",
  -macro_pars$Sigma[2, 2] / (sqrt(macro_pars$Alpha[2, 2]) * constraints$beta),
  "$H_2$",
  "$V^*_{g_\\mathrm{s,max}}$",
  macro_pars$Sigma[2, 2] / (2 * sqrt(macro_pars$Alpha[2, 2])),
  "$H_3$",
  "$V^*_{d_\\mathrm{S}}$",
  constraints$beta^2 * Vstar_h3[1, 1] + Vstar_h3[2, 2],
  "$H_3$",
  "$V^*_{a_\\mathrm{S}}$",
  Vstar_h3[1, 1],
  "$H_3$",
  "$V^*_{d_\\mathrm{S},a_\\mathrm{S}}$",
  -Vstar_h3[1, 2] - constraints$beta * Vstar_h3[1, 1],
  "$H_3$",
  "$V^*_{g_\\mathrm{s,max}}$",
  Vstar_h3[2, 2]
) |>
  mutate(component1 = factor(component1, levels = components))

# Simulation results ----
c("df_Mratio", "df_rM", "df_omega", "df_phi") |>
  walk(\(.x) {
    assign(
      .x,
      read_rds(glue("objects/{.x}.rds")) |>
        mutate(
          hypothesis1 = case_when(hypothesis == "h2" ~ "$H_2$", hypothesis == "h3" ~ "$H_3$", ),
          component1 = factor(
            case_when(
              component == "V_aS" ~ "$V^*_{a_\\mathrm{S}}$",
              component == "V_dS" ~ "$V^*_{d_\\mathrm{S}}$",
              component == "V_dSaS" ~ "$V^*_{d_\\mathrm{S},a_\\mathrm{S}}$",
              component == "V_gmax" ~ "$V^*_{g_\\mathrm{s,max}}$"
            ),
            levels = components
          )
        ) |>
        pivot_wider(names_from = quantity, values_from = value),
      envir = .GlobalEnv
    )
  })

# 1. mRatio ----

gp_mRatio = ggplot(df_Mratio, aes(Mratio, estimate, color = component1)) +
  facet_grid(hypothesis1 ~ ., scales = "free_y") +
  geom_hline(
    data = df_eq,
    aes(yintercept = estimate, color = component1),
    linetype = 2,
    linewidth = 1.1,
    alpha = 0.5
  ) +
  geom_pointinterval(mapping = aes(ymin = lower, ymax = upper), alpha = 0.5) +
  scale_x_log10() +
  scale_color_discrete(name = NULL) +
  xlab("Mutational ratio ($M_{d_\\mathrm{S}} / M_{a_\\mathrm{S}}$)") +
  ylab("Stationary (co)variance")

# 2. rho_M ----
gp_rM = ggplot(df_rM, aes(rho_M, estimate, color = component1)) +
  facet_grid(hypothesis1 ~ ., scales = "free_y") +
  geom_hline(
    data = df_eq,
    aes(yintercept = estimate, color = component1),
    linetype = 2,
    linewidth = 1.1,
    alpha = 0.5
  ) +
  geom_pointinterval(mapping = aes(ymin = lower, ymax = upper), alpha = 0.5) +
  scale_color_discrete(name = NULL) +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5)) +
  xlab("Mutational correlation ($\\rho_\\mathrm{M}$)") +
  ylab("Stationary (co)variance")

# 3. omega ----
gp_omega = ggplot(df_omega, aes(omega, estimate, color = component1)) +
  facet_grid(hypothesis1 ~ ., scales = "free_y") +
  geom_hline(
    data = df_eq,
    aes(yintercept = estimate, color = component1),
    linetype = 2,
    linewidth = 1.1,
    alpha = 0.5
  ) +
  geom_pointinterval(mapping = aes(ymin = lower, ymax = upper), alpha = 0.5) +
  scale_x_log10() +
  scale_color_discrete(name = NULL) +
  xlab("Strength of selection on $g_\\mathrm{s,max}$ ($\\omega$)") +
  ylab("Stationary (co)variance")

# 4. phi ----
gp_phif = ggplot(df_phi, aes(phi_f, estimate, color = component1)) +
  facet_grid(hypothesis1 ~ ., scales = "free_y") +
  geom_hline(
    data = df_eq,
    aes(yintercept = estimate, color = component1),
    linetype = 2,
    linewidth = 1.1,
    alpha = 0.5
  ) +
  geom_pointinterval(mapping = aes(ymin = lower, ymax = upper), alpha = 0.5) +
  scale_x_log10(breaks = c(0.5, 1, 2)) +
  scale_color_discrete(name = NULL) +
  xlab("Strength of selection on $f_\\mathrm{S}$ ($\\phi_f$)") +
  ylab("Stationary (co)variance")

gp_phia = ggplot(filter(df_phi, hypothesis == "h3"),
                 aes(phi_a, estimate, color = component1)) +
  geom_hline(
    data = filter(df_eq, hypothesis1 == "$H_3$"),
    aes(yintercept = estimate, color = component1),
    linetype = 2,
    linewidth = 1.1,
    alpha = 0.5
  ) +
  geom_pointinterval(mapping = aes(ymin = lower, ymax = upper), alpha = 0.5) +
  scale_x_log10(breaks = c(0.5, 1, 2)) +
  scale_color_discrete(name = NULL) +
  xlab("Strength of selection on $a_\\mathrm{S}$ ($\\phi_a$)") +
  ylab("Stationary (co)variance")

# Write figures ----

tikz(file = "ms/figures/fig_mRatio.tex", width = 4, height = 4)
print(gp_mRatio)  
dev.off()

tikz(file = "ms/figures/fig_rM.tex", width = 4, height = 4)
print(gp_rM)  
dev.off()

tikz(file = "ms/figures/fig_omega.tex", width = 4, height = 4)
print(gp_omega)  
dev.off()

tikz(file = "ms/figures/fig_phif.tex", width = 4, height = 4)
print(gp_phif)  
dev.off()

tikz(file = "ms/figures/fig_phia.tex", width = 4, height = 3)
print(gp_phia)  
dev.off()
