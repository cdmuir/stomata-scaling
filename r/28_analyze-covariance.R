# Code to analyze range of possible covariance between stomatal density and size when fs is bounded.
source("r/header.R")

dat = read_rds("objects/dat.rds") |>
  # convert d_s from log(mm^-2) to log(um^-2)
  mutate(d_s = d_s - log(1e6))

mu_ds = mean(dat$d_s)
mu_as = mean(dat$a_s)
sd_ds = sd(dat$d_s)
sd_as = sd(dat$a_s)

# Variables
fSmax = tibble(name = c("1", "1/2", "1/3"), value = c(1, 0.5, 1/3), z_max = log(value))

variable_list = crossing(
  rho = seq(-0.9, 0.9, length.out = 100),
  z_max = fSmax$z_max
) |>
  as.list()

plan(multisession, workers = 10)

df_rho = future_map2_dfr(
  variable_list$rho,
  variable_list$z_max,
    \(rho, z_max, Mu, sd_ds, sd_as, B) {
      
      # Define covariance matrix
      Sigma = matrix(c(sd_ds ^ 2, rho * sd_ds * sd_as, rho * sd_ds * sd_as, sd_as ^ 2), nrow = 2)
      
      # numerical integration
      num = get_cov(Mu, Sigma, B, z_max)
      
      # simulation
      sim = mvnfast::rmvn(1e7, Mu, Sigma) |>
        as_tibble(.name_repair = "unique") |>
        set_colnames(c("d_s", "a_s")) |>
        mutate(z = B[1] * d_s + B[2] * a_s) |>
        filter(z <= z_max) |>
        summarize(
          mu_ds = mean(d_s),
          mu_as = mean(a_s),
          var_ds = var(d_s),
          var_as = var(a_s),
          cov = cov(d_s, a_s)
        )
      
      tibble(
        rho_actual = rho,
        z_max = z_max,
        type = c("actual", "numerical", "simulated"),
        cov = c(rho * (sd_ds * sd_as), num["cov"], sim$cov),
        sd_ds = c(sd_ds, sqrt(num["var1"]), sqrt(sim$var_ds)),
        sd_as = c(sd_as, sqrt(num["var2"]), sqrt(sim$var_as)),
        mu_ds = c(Mu[1], num["mu1"], sim$mu_ds),
        mu_as = c(Mu[2], num["mu2"], sim$mu_as)
      )
      
    },
    Mu = c(mu_ds, mu_as),
    sd_ds = sd_ds,
    sd_as = sd_as,
    B = c(1, 1),
  .progress = TRUE
  )

df_rho |> 
  left_join(select(fSmax, name, z_max), by = join_by(z_max)) |>
  mutate(rho_apparent = cov / (sd_ds * sd_as)) |>
  filter(type != "actual") |>
  ggplot(aes(rho_actual, rho_apparent, color = name)) +
  facet_grid(~type) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "True correlation", y = "Apparent correlation") +
  scale_color_discrete(name = expression(italic(f)[paste(S,",",max)])) +
  theme_cowplot()

ggsave("figures/20_analyze-covariance.png", width = 8, height = 4, dpi = 300)