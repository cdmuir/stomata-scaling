# It would be nice to calculate 'scale' so that ellipse was exactly adjacent to
# empirical boundary. Currently, I tweaked values visually.
#
# Also need to label geometric constraint lines in figure
#
source("r/header.R")
dat = read_rds("objects/dat.rds")

mu_ds = mean(dat$d_s)
mu_as = mean(dat$a_s)
var_ds = var(dat$d_s)
var_as = var(dat$a_s)

df_constraint = crossing(f_s = c(exp(max(dat$log_fs)), 1), A_s = exp(seq(min(dat$a_s), max(dat$a_s), length.out = 101))) |>
  mutate(
    D_s = f_s * 1e6 / A_s,
    a_s = log(A_s),
    d_s = log(D_s),
    f_s = set_units(D_s, 1 / mm^2) * set_units(A_s, um^2)
  )

# Panel A -----
# Many relationships between stomatal size and density are geometrically possible

df_ellipse = list(
  correlation = c("negative", "zero", "positive"),
  rho = c(-0.5, 0, 0.5),
  scale = c(0.75, 0.55, 0.45)
) |>
  pmap_dfr( ~ {
    ellipse(
      x = matrix(c(1, ..2, ..2, 1), ncol = 2),
      scale = c(..3, ..3),
      centre = c(mu_as, mu_ds),
      level = 0.95
    ) |>
      as_tibble(.name_repair = NULL) |>
      rename(a_s = x, d_s = y) |>
      mutate(correlation = ..1,
             rho = ..2,
             scale = ..3)
    
  }) |>
  mutate(correlation = factor(correlation, levels = c("negative", "zero", "positive")))

panel_a = ggplot(df_ellipse, aes(exp(a_s), exp(d_s))) +
  facet_grid(. ~ correlation) +
  scale_x_continuous(
    trans = "log",
    limits = exp(mu_as + c(-4.5, 4.5) * sqrt(var_as)),
    breaks = c(10, 100, 1000)
  ) +
  scale_y_continuous(
    trans = "log",
    limits = exp(mu_ds + c(-4.5, 5) * sqrt(var_as)),
    breaks = c(10, 100, 1000)
  ) +
  geom_polygon(fill = "grey") +
  geom_line(data = df_constraint, mapping = aes(linetype = as.factor(f_s)), ) +
  xlab(expression(paste("Stomatal size [µ", m^2, "], log-scale"))) +
  ylab(expression(atop(
    Stomatal ~ density, paste(group("[", pores ~ mm^-2, "]"), ",", ~ log -
                                scale)
  ))) +
  ggtitle("Many relationships between stomatal size and\ndensity are geometrically possible") +
  coord_equal() +
  theme_cowplot() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "none",
    plot.title.position = "plot",
    strip.text = element_text(size = 12)
  ) +
  NULL

# Panel B ----
assumptions = c(
  "atop(Stabilizing~selection,on~optimal~italic(g)[paste(s,',',max)])",
  "atop(Directional~selection,to~minimize~italic(f)[S])",
  "atop(Stabilizing~selection~on,optimal~stomatal~size)"
)

df_area = crossing(x = seq(-4.5, 4.5, 0.1), nesting(
  assumption = factor(assumptions, levels = assumptions),
  mu = c(0, -4.5, 0)
)) |>
  mutate(y = dnorm(x, mean = mu, sd = 1))

df_arrow = bind_rows(
  tibble(
    assumption = assumptions[1],
    group = letters[seq_len(2)],
    x = c(-4.5, 4.5),
    y = rep(0.2, 2),
    xend = c(-1.5, 1.5),
    yend = y
  ),
  tibble(
    assumption = assumptions[2],
    group = letters[seq_len(2)],
    x = 0,
    y = 0.2,
    xend = -3,
    yend = y
  ),
  tibble(
    assumption = assumptions[3],
    group = letters[seq_len(2)],
    x = c(-4.5, 4.5),
    y = rep(0.2, 2),
    xend = c(-1.5, 1.5),
    yend = y
  )
) |>
  mutate(assumption = factor(assumption, levels = assumptions))


panel_b = ggplot(df_area, aes(x, y)) +
  facet_grid(. ~ assumption, labeller = label_parsed) +
  geom_area(fill = "darkgrey") +
  geom_segment(
    data = df_arrow,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      group = group
    ),
    arrow = arrow(length = unit(0.1, "inches"))
  ) +
  xlab("Trait value") +
  ylab("Fitness") +
  ggtitle("Assumptions underlying competing hypotheses for\ninverse size-density scaling") +
  theme_cowplot() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 12),
    legend.position = "none",
    plot.title.position = "plot",
    strip.text = element_text(size = 12)
  ) +
  NULL

plot_grid(
  panel_a,
  panel_b,
  ncol = 1,
  align = "v",
  axis = "l",
  rel_heights = c(1, 1),
  labels = "auto"
)

ggsave(
  "ms/figures/hypotheses.png",
  width = 6.5,
  height = 6.5,
  units = "in"
)
