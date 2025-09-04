source("r/header.R")

fit_phylolm1 = read_rds("objects/fit_phylolm1.rds")
fit_phylolm2 = read_rds("objects/fit_phylolm2.rds")
dat = read_rds("objects/dat.rds")

mu_ds = mean(dat$d_s)
mu_as = mean(dat$a_s)
mu_fs = mean(dat$log_fs)
mu_gmax = mean(dat$log_gmax)
var_ds = var(dat$d_s)
var_as = var(dat$a_s)
var_fs = var(dat$log_fs)
var_gmax = var(dat$log_gmax)

boot = fit_phylolm2$bootstrap |>
  as.data.frame() %>%
  mutate(sim = 1:nrow(.)) |>
  rename(b_as = a_s)

s_ranges = dat |>
  filter(species %in% rownames(fit_phylolm2$X)) |>
  group_by(group) |>
  summarize(
    min_as = min(a_s),
    max_as = max(a_s),
    med_as = mean(a_s),
    med_ds = mean(d_s),
    .groups = "drop"
  ) |>
  ungroup()

df_plot = crossing(
  a_s = seq(min(s_ranges$min_as), max(s_ranges$max_as), 
               length.out = 1e2),
  sim = 1:nrow(boot)
) |>
  full_join(boot, by = "sim") |>
  mutate(
    y_Angiosperms = `(Intercept)` + b_as * a_s,
    y_Gymnosperms = `(Intercept)` + groupGymnosperms +
      c(b_as + `a_s:groupGymnosperms`) * a_s,
    y_Pteridophytes = `(Intercept)` + groupPteridophytes +
      c(b_as + `a_s:groupPteridophytes`) * a_s
  ) |>
  group_by(a_s) |>
  point_interval(
    y_Angiosperms, y_Gymnosperms, y_Pteridophytes,
    .width = 0.95, .point = median, .interval = qi
  ) |>
  ungroup() |>
  gather(key, y, -a_s, -.interval, -.point, -.width) |>
  mutate(
    group = str_replace(key, "^y_([[:alpha:]]+).*[[:alpha:]]*$", "\\1"),
    quantity = str_replace(key, "^y_[[:alpha:]]+.([[:alpha:]]*)$", "\\1"),
    quantity = ifelse(quantity == "", "estimate", quantity)
    ) |>
  full_join(select(s_ranges, group, min_as, max_as), by = "group") |>
  filter(
    a_s >= min_as,
    a_s <= max_as
  ) |>
  select(a_s, y, group, quantity) |>
  spread(quantity, y)

df_pred = df_plot |>
  select(a_s, group) |>
  crossing(slope = c(-0.5, -1)) |>
  full_join(select(s_ranges, group, med_as, med_ds), by = "group") |>
  mutate(intercept = med_ds - slope * med_as) |>
  mutate(estimate = slope * a_s + intercept)

df_text = boot |>
  transmute(
    slope_angiosperms = b_as,
    slope_gymnosperms = b_as + `a_s:groupGymnosperms`,
    slope_pteridophytes = b_as + `a_s:groupPteridophytes`
  ) |>
  gather(parameter, value) |>
  mutate(value = -1 * value) |>
  group_by(parameter) |>
  point_interval(.width = 0.95, .point = median, .interval = qi) |>
  select(parameter, value, .lower, .upper) |>
  mutate(sl = ifelse(.lower < 0, "-", ""), su = ifelse(.upper < 0, "-", "")) |>
  mutate_if(is.numeric, abs) |>
  mutate_if(is.numeric, formatC, digits = 2, flag = "#") |>
  mutate_if(is.numeric, as.character) |>
  mutate(
    text = glue("hat(beta)==`{estimate}`~({sl}`{lower}`~`–`~{su}`{upper}`)",
                estimate = value, sl = sl, lower = .lower, 
                su = su, upper = .upper),
    group = parameter |>
      str_replace("^slope_([[:alpha:]]+)$", "\\1") |>
      str_to_sentence(),
    x = exp(mu_as - 4.5 * sqrt(var_as)),
    y = exp(mu_ds + 5 * sqrt(var_as))
  )

gp1 = ggplot(dat, aes(exp(a_s), exp(d_s))) +
  scale_x_continuous(
    trans = "log", 
    limits = exp(mu_as + c(-4.5, 4.5) * sqrt(var_as)),
    # limits = c(10, 3100), 
    breaks = c(10, 100, 1000)
  ) + 
  scale_y_continuous(
    trans = "log", 
    limits = exp(mu_ds + c(-4.5, 5) * sqrt(var_as)),
    # limits = c(10, 2100), 
    breaks = c(10, 100, 1000)
  ) + 
  facet_grid(~ group) +
  geom_point(
    data = dplyr::select(dat, -group), 
    color = "grey50"
  ) +
  geom_point(alpha = 0.5, size = 2, shape = 19) +
  # geom_ribbon(
  #   data = rename(df_plot, d_s = estimate), 
  #   mapping = aes(ymin = exp(lower), ymax = exp(upper)),
  #   fill = "black", color = NA, alpha = 0.5#, inherit.aes = FALSE
  # ) +
  # geom_line(
  #   data = df_pred, 
  #   mapping = aes(y = exp(estimate), linetype = as.factor(slope)),
  #   color = "grey", linewidth = 1.1, lineend = "round"
  # ) +
  # geom_line(
  #   data = df_plot, 
  #   mapping = aes(y = exp(estimate)),
  #   color = "tomato", size = 1.5, lineend = "round"
  # ) +
  # geom_text(
  #   data = df_text,
  #   mapping = aes(x, y, label = text),
  #   parse = TRUE,
  #   hjust = 0,
  #   vjust = 1
  # ) +
  xlab(expression(paste("Stomatal size [μ", m ^2, "], log-scale"))) +
  ylab(expression(atop(Stomatal~density, paste(group("[", pores~mm ^-2, "]"), ",",~log-scale)))) +
  coord_equal() +
  theme_bw() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 12),
    legend.position = "none"
  ) +
  NULL

gp2 = ggplot(dat, aes(exp(a_s), exp(log_gmax))) +
  scale_x_continuous(
    trans = "log", 
    limits = exp(mu_as + c(-4.5, 4.5) * sqrt(var_as)),
    # limits = c(10, 3100), 
    breaks = c(10, 100, 1000)
  ) + 
  scale_y_continuous(
    trans = "log", 
    limits = exp(mu_gmax + c(-4.5, 5) * sqrt(var_as)),
    # limits = c(10, 2100), 
    breaks = c(0.1, 1, 10)
  ) + 
  facet_grid(~ group) +
  geom_point(
    data = dplyr::select(dat, -group), 
    color = "grey50"
  ) +
  geom_point(alpha = 0.5, size = 2, shape = 19) +
  xlab(expression(paste("Stomatal size [μ", m ^2, "], log-scale"))) +
  ylab(expression(atop(maximum~stomatal,conductance~(italic(g)[paste(s,,",",max)])))) +
  coord_equal() +
  theme_bw() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 12),
    legend.position = "none"
  ) +
  NULL

gp3 = ggplot(dat, aes(exp(d_s), exp(log_gmax))) +
  scale_x_continuous(
    trans = "log", 
    limits = exp(mu_ds + c(-4.5, 5) * sqrt(var_as)),
    # limits = c(10, 3100), 
    breaks = c(10, 100, 1000)
  ) + 
  scale_y_continuous(
    trans = "log", 
    limits = exp(mu_gmax + c(-4.5, 5) * sqrt(var_as)),
    # limits = c(10, 2100), 
    breaks = c(0.1, 1, 10)
  ) + 
  facet_grid(~ group) +
  geom_point(
    data = dplyr::select(dat, -group), 
    color = "grey50"
  ) +
  geom_point(alpha = 0.5, size = 2, shape = 19) +
  xlab(expression(atop(Stomatal~density, paste(group("[", pores~mm ^-2, "]"), ",",~log-scale)))) +
  ylab(expression(atop(maximum~stomatal,conductance~(italic(g)[paste(s,",",max)])))) +
  coord_equal() +
  theme_bw() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 12),
    legend.position = "none"
  ) +
  NULL

gp4 = ggplot(dat, aes(exp(log_fs), exp(log_gmax))) +
  scale_x_continuous(
    trans = "log",
    limits = exp(mu_fs + c(-5.5, 4.5) * sqrt(var_fs)),
    breaks = c(0.005, 0.05, 0.5)
  ) +
  scale_y_continuous(
    trans = "log", 
    limits = exp(mu_gmax + c(-5, 5) * sqrt(var_gmax)),
    breaks = c(0.1, 1, 10)
  ) + 
  facet_grid(~ group) +
  geom_point(
    data = dplyr::select(dat, -group), 
    color = "grey50"
  ) +
  geom_point(alpha = 0.5, size = 2, shape = 19) +
  # xlab("fraction of epidermal area allocated to stomata ($f_\\mathrm{S}$)") +
  # ylab("maximum stomatal conductance ($g_\\mathrm{max}$)") +
  xlab(expression(fraction~of~epidermal~area~allocated~to~stomata~(italic(f)[S]))) +
  ylab(expression(atop(maximum~stomatal,conductance~(italic(g)[paste(s,',',max)])))) +
  coord_equal() +
  theme_bw() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 12),
    legend.position = "none"
  ) +
  NULL

gp5 = plot_grid(gp1, gp3, gp2, ncol = 1, labels = "auto")

ggsave(
  "ms/figures/results.png",
  plot = gp5,
  width = 6.5,
  height = 9,
  units = "in"
)

ggsave(
  file = "ms/figures/fig-gmax-fs.pdf",
  plot = gp4,
  width = 6.5,
  height = 3.25,
  units = "in"
)
