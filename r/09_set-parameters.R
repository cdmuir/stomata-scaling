# Set parameters for simulations
source("r/header.R")

dat = read_rds("objects/dat.rds")
phy = read_rds("objects/phy.rds")
b = biophysical_constant(2.49e-5, 2.24e-2)
m = morphological_constant(0.5, 0.5, 0.5) 
theta_gmax = mean(dat$log_gmax)

# Fit to get macroevolutionary parameters for g_max, d_S, and a_S
fit_gmax = phylolm(
  formula = log_gmax ~ group + grass + shrub + tree,
  data = column_to_rownames(dat, "species"),
  phy = phy,
  model = "OUrandomRoot",
  upper.bound = 10
)

fit_dS = phylolm(
  d_s ~ group + grass + shrub + tree,
  column_to_rownames(dat, "species"),
  phy,
  model = "OUrandomRoot",
  upper.bound = 10
)

fit_aS = phylolm(
  a_s ~ group + grass + shrub + tree,
  column_to_rownames(dat, "species"),
  phy,
  model = "OUrandomRoot",
  upper.bound = 10
)

z0_fS = log(1e-6) # this is for unit conversion
z0_gmax = log(b * m)

# Hyperparameters

## Minimum fS for h2
fS_min = log(1e-3)

## Scaling parameter
## Assume that beta = 0.5 for simulations
beta = 0.5

## Maximum absolute fitness 
## this is not important because selection is based on relative fitness
W_max = 1

## Stationary variance in g_max
V_gmax_stat = fit_gmax$sigma2 / (2 * fit_gmax$optpar)

## Proportion of environmental variance relative to stationary variance in Z
p_Ve = 0.01

## Proportion of mutation variance relative to environmental variance
p_Vm = 0.001

## Proportion of selection variance in g_max relative to environmental variance
p_omega_gmax = 10

## Proportion of selection variance in f_S relative to selection variance in g_max
phi_f = 1

## Proportion of selection variance in a_S relative to selection variance in g_max
phi_a = 1

## Mutational correlation between log(density) and log(size) traits
r_da = 0

## Macroevolutionary variances
# Convert units from my to years (assume 1 generation per year)
# a_S ("stomatal size")
theta_aS = coef(fit_aS)["(Intercept)"]
alpha_aS = fit_aS$optpar / 1e6
sigma_aS = sqrt(fit_aS$sigma2 / 1e6)
# g_max
theta_gmax = coef(fit_gmax)["(Intercept)"] 
alpha_gmax = fit_gmax$optpar / 1e6
sigma_gmax = sqrt(fit_gmax$sigma2 / 1e6)

## Macroevolutionary correlations
r_alpha_aS_gmax = 0
r_sigma_aS_gmax = 0

# Microevolutionary parameters ----
micro_pars = list(
  V_e = p_Ve * V_gmax_stat # was 5e-3, # 1% of Stationary variance in Z
)

# 5e-6, # V_m = V_e / 1000
micro_pars$M = matrix(c(
  p_Vm * micro_pars$V_e,
  rep(r_da * p_Vm * micro_pars$V_e, 2),
  p_Vm * micro_pars$V_e
), 2, 2)
micro_pars$omega = p_omega_gmax * micro_pars$V_e # 5e-2 # V_s = 10 V_e
micro_pars$phi_a = phi_a 
micro_pars$phi_f = phi_f
micro_pars$r_da = r_da

# Macroevolutionary parameters ----
macro_pars = list(
  Theta = matrix(c(theta_aS, theta_gmax), 2, 1),
  Alpha = matrix(c(alpha_aS ^ 2, rep(r_alpha_aS_gmax * alpha_aS * alpha_gmax, 2), alpha_gmax ^ 2), 2, 2),
  Sigma = matrix(c(sigma_aS ^ 2, rep(r_sigma_aS_gmax * sigma_aS * sigma_gmax, 2), sigma_gmax ^ 2), 2, 2)
)

# check that the covariance matrix is specified realistically
cov_MVOU1(1e7, sqrtm(macro_pars$Alpha), macro_pars$Sigma)

# Constraints ----

constraints = list(
  z0_fS = z0_fS,
  z0_gmax = z0_gmax,
  fS_min = fS_min,
  beta = beta
)

## List of hyperparameters to export
write_rds(list(
  W_max = W_max,
  V_gmax_stat = V_gmax_stat,
  p_Ve = p_Ve,
  p_Vm = p_Vm,
  phi_a = phi_a,
  phi_f = phi_f,
  p_omega_gmax = p_omega_gmax,
  r_da = r_da,
  theta_aS = theta_aS,
  alpha_aS = alpha_aS,
  sigma_aS = sigma_aS,
  theta_gmax = theta_gmax,
  alpha_gmax = alpha_gmax,
  sigma_gmax = sigma_gmax,
  r_sigma_aS_gmax = r_sigma_aS_gmax,
  r_sigma_aS_gmax = r_sigma_aS_gmax
), "objects/hyperparameters.rds")

## List of parameters to export

write_rds(
  list(
    constraints = constraints,
    micro_pars = micro_pars,
    macro_pars = macro_pars
  ),
  "objects/parameters.rds"
)
