# Predicted interspecific OLS and SMA slopes from macro parameters
slope_h34 = function(macro_pars, constraints) {
  
  V = StationaryVariance(macro_pars$Alpha, macro_pars$Sigma)
  
  V_aS = V[1,1]
  V_gmax = V[2,2]
  Cov_aS_gmax = V[1,2]
  V_dS = V_gmax + constraints$beta ^ 2 * V_aS - 2 * constraints$beta * Cov_aS_gmax
  Cov_dS_aS = -constraints$beta * V_aS
  
  list(
    ols = Cov_dS_aS / V_dS,
    sma = -sqrt(V_aS / V_dS)
  )
}

get_del_logW = function(..., model) {
  switch(
    model,
    "h1" = get_del_logW_h1(...),
    "h2" = get_del_logW_h2(...),
    "h3" = get_del_logW_h3(...)
  )
}

# For H1, see solve-h1.py
get_del_logW_h1 = function(g_opt, da_bar, omega, z0_gmax, beta) {
  
  dS = da_bar[1]
  aS = da_bar[2]
  
  c(
   -(z0_gmax + dS + beta * aS - g_opt)/ omega,
   -beta * (z0_gmax + dS + beta * aS - g_opt) / omega
  )
  
}

# For H2, see solve-h2.py
get_del_logW_h2 = function(g_opt, fS_min, da_bar, phi_f, omega, z0_fS, z0_gmax, beta) {
  
  dS = da_bar[1]
  aS = da_bar[2]

  c(
    # d log(W) / d dS_bar
    -(z0_gmax + dS + beta * aS - g_opt) / omega - (z0_fS + aS + dS - fS_min) / (omega * phi_f),
    # d log(W) / d aS_bar
    -beta * (z0_gmax + dS + beta * aS - g_opt) / omega - (z0_fS + aS + dS - fS_min) / (omega * phi_f)  
  )
 
}

# For H3, see solve-h3.py
get_del_logW_h3 = function(ag_opt, da_bar, phi_a, phi_f, omega, z0_gmax, beta) {
  
  dS = da_bar[1]
  aS = da_bar[2]
  aS_opt = ag_opt[1]
  g_opt = ag_opt[2]
  
  c(
    # d log(W) / d dS_bar
    (-phi_a * phi_f * (2 * aS * beta + 2 * dS - 2 * g_opt + 2 * z0_gmax) - phi_a) /
      (2 * omega * phi_a * phi_f),
    # d log(W) / d aS_bar
    (
      -2 * beta * phi_a * phi_f * (aS * beta + dS - g_opt + z0_gmax) - phi_a - phi_f *
        (2 * aS - 2 * aS_opt)
    ) / (2 * omega * phi_a * phi_f)
  )
  
}

# Approximate change in mean dS and aS in one generation
approx_dda_bar = function(ag_opt, da_bar, micro_pars, constraints, G, model) {
  switch(
    model,
    "h1" = approx_dda_bar_h1(ag_opt, da_bar, micro_pars, constraints, G),
    "h2" = approx_dda_bar_h2(ag_opt, da_bar, micro_pars, constraints, G),
    "h3" = approx_dda_bar_h3(ag_opt, da_bar, micro_pars, constraints, G)
  )
}

approx_dda_bar_h1 = function(ag_opt, da_bar, micro_pars, constraints, G) {
  
  del_logW = get_del_logW(
    g_opt = ag_opt[2],
    da_bar = da_bar,
    omega = micro_pars$omega,
    z0_gmax = constraints$z0_gmax,
    beta = constraints$beta,
    model = "h1"
  )
  
  G %*% del_logW
  
}

approx_dda_bar_h2 = function(ag_opt, da_bar, micro_pars, constraints, G) {
  
  del_logW = get_del_logW(
    g_opt = ag_opt[2],
    fS_min = constraints$fS_min,
    da_bar = da_bar,
    phi_f = micro_pars$phi_f,
    omega = micro_pars$omega,
    z0_fS = constraints$z0_fS,
    z0_gmax = constraints$z0_gmax,
    beta = constraints$beta,
    model = "h2"
  )
  
  G %*% del_logW
  
}

approx_dda_bar_h3 = function(ag_opt, da_bar, micro_pars, constraints, G) {
  
  del_logW = get_del_logW(
    ag_opt = ag_opt,
    da_bar = da_bar,
    phi_a = micro_pars$phi_a,
    phi_f = micro_pars$phi_f,
    omega = micro_pars$omega,
    z0_gmax = constraints$z0_gmax,
    beta = constraints$beta,
    model = "h3"
  )
  
  G %*% del_logW
  
}

# Set up output
setup_output = function(ag_opt, n_gen, thin, type) {
  
  # 1. Set up output
  i = seq(1, n_gen, by = thin)
  tibble(
    gen = i,
    # trait optima
    aS_opt = ag_opt[i,1],
    g_opt = ag_opt[i,2],
    # mean trait values 
    dag_bar = vector("list", length(i)),
    ddag_bar = vector("list", length(i)),
    # additive genetic covariance
    G = vector("list", length(i)),
    dG = vector("list", length(i)),
    type = type
  )
  
}

# Simulate selection for one lineage for n_gen generations using approximate equations for change in mean trait values
selection_approximate = function(n_gen,
                                 ag_opt,
                                 init,
                                 constraints,
                                 micro_pars,
                                 macro_pars,
                                 model,
                                 drift,
                                 thin,
                                 seed) {
  set.seed(seed)
  
  # 1. Set up output
  out = setup_output(ag_opt, n_gen, thin, type = "approximate") 
  
  dag_bar = init$dag_bar
  G = init$G
  
  # Curvature and orientation of fitness surface
  H = get_H(micro_pars$phi_a, micro_pars$omega, constraints$beta, model = model)
  
  # 2. Selection, mutation, and drift
  for (t in seq_len(n_gen)) {
    
    g_bar = constraints$z0_gmax + dag_bar[1,1] + constraints$beta * dag_bar[2,1]
    dda_bar = approx_dda_bar(
      ag_opt[t,],
      dag_bar[1:2,1],
      micro_pars,
      constraints,
      G,
      model = model
    )
    
    # Gaussian approximation for change in additive genetic variance
    G1 = G + G %*% H %*% G
    
    if (drift) {
      dda_bar = dda_bar + t(mvnfast::rmvn(1, rep(0, 2), sigma = G / micro_pars$N))
      G1 = rWishart(n_sim, df = micro_pars$N - 1, Sigma = G1) / micro_pars$N
    }
    
    G1 = G1 + micro_pars$M 
    dg_bar = dda_bar[1,1] + constraints$beta * dda_bar[2,1]
    ddag_bar = matrix(c(dda_bar, dg_bar), 3, 1)
    
    # 3. Record statistics 
    if (t %in% out$gen) {
      ii = which(out$gen == t)
      out$dag_bar[[ii]]  = dag_bar
      out$ddag_bar[[ii]] = ddag_bar
      out$G[[ii]]  = G
      out$dG[[ii]] = G1 - G
    }
    
    # 4. Update
    dag_bar = dag_bar + ddag_bar
    G = G1
    
  }
  
  out
  
}

# Simulate lineage with moving optima and selection on a composite trait
# 
# n_gen (integer(1)): Number of generations to simulate
# constraints (list): Parameters setting constraints on evolutionary dynamics
# micro_pars (list): Parameters for the microevolutionary model
# macro_pars (list): Parameters for the macroevolutionary model
# model (character(1)): Which model to use for the microevolutionary model? Options are "h1", "h2", and "h3"
# individual_based (flag): Should individual based simulations be performed?
# drift (flag): Should drift be added to the approximate model?
# thin (integer(1)): Thinning parameter for output
# init_method (character(1)): How should the initial z_bar be set?
# custom_init (list): Custom initial conditions for the simulation (not checked)
# n_ind_sim (integer(1)): Number of individual-based simulations
# S_init (matrix(2, 2)): Initial covariance matrix for the macroevolutionary model. Only needed if init_method == "random".
simulate_lineage = function(n_gen,
                            constraints,
                            micro_pars,
                            macro_pars,
                            x = NULL,
                            model,
                            individual_based,
                            drift,
                            thin,
                            init_method,
                            custom_init = NULL,
                            n_ind_sim = 1,
                            S_init = NULL,
                            seed = NULL) {
  
  if (is.null(custom_init)) {
    if (is.null(x)) {
      stop("If custom_init is NULL, x must be provided. x is initial guess for G* in nleqslv.")
    }
    init = initialize_sim(n_gen,
                          constraints,
                          micro_pars,
                          macro_pars,
                          x,
                          model,
                          init_method,
                          S_init,
                          seed)
  } else {
    init = custom_init
  }
  
  out = selection_approximate(n_gen,
                              init$ag_opt,
                              init$micro_init,
                              constraints,
                              micro_pars,
                              macro_pars,
                              model,
                              drift,
                              thin,
                              seed) |>
    mutate(sim = 1)
  if (individual_based) {
    for (s in seq_len(n_ind_sim)) {
      init$micro_init$xy = mvnfast::rmvn(micro_pars$N,
                                         init$micro_init$xyz_bar[1:2, 1],
                                         sigma = init$micro_init$Sigma_Am)
      
      out = bind_rows(
        out,
        selection_individual(
          n_gen,
          init$yz_opt,
          init$micro_init,
          constraints,
          micro_pars,
          macro_pars,
          thin,
          seed
        ) |>
          mutate(sim = s)
      )
    }
  }
  
  out |>
    mutate(across(where(is.list), ~ map(.x, as.vector))) |>
    unnest_wider(where(is.list), names_sep = "_")
  
}

initialize_sim = function(n_gen,
                          constraints,
                          micro_pars,
                          macro_pars,
                          x = NULL,
                          model,
                          init_method,
                          S_init = NULL,
                          seed) {
  
  set.seed(seed)
  init_method = switch(init_method,
                       random = "random",
                       stationary = "stationary",
                       Theta = "Theta")
  
  # 1. Initial mean trait values
  dag_bar0 = get_initial_mean_trait(macro_pars, constraints, init_method, S_init)
  
  # 2. Random trajectory of trait optima
  ag_opt = rTrajectoryMVOU(n_gen + 1,
                           dag_bar0[2:3, 1],
                           1,
                           macro_pars$Alpha,
                           macro_pars$Theta[,1],
                           macro_pars$Sigma,
                           seed)[2:(n_gen + 1),]
  
  # 3. Start population genetic covariance at mutation-selection balance
  micro_init = get_micro_init(dag_bar0, micro_pars, constraints, x = x, model = model)
  
  list(ag_opt = ag_opt, micro_init = micro_init)
  
}

get_initial_mean_trait = function(macro_pars, constraints, init_method, S_init = NULL) {
  
  if (init_method == "random") {
    ag_bar0 = rmvn(1, macro_pars$Theta, S_init)
    dag_bar0 = matrix(c(
      ag_bar0[1, 2] - constraints$z0_gmax - constraints$beta * ag_bar0[1, 1],
      ag_bar0[1, 1],
      ag_bar0[1, 2]
    ),
    3,
    1)
  }
  if (init_method == "stationary") {
    S = StationaryVariance(macro_pars$Alpha, macro_pars$Sigma)
    ag_bar0 = rmvn(1, macro_pars$Theta, S)
    dag_bar0 = matrix(c(
      ag_bar0[1, 2] - constraints$z0_gmax - constraints$beta * ag_bar0[1, 1],
      ag_bar0[1, 1],
      ag_bar0[1, 2]
    ),
    3,
    1)
  }
  if (init_method == "Theta") {
    dag_bar0 = matrix(
      c(
        macro_pars$Theta[2, 1] - constraints$z0_gmax - constraints$beta * macro_pars$Theta[1, 1],
        macro_pars$Theta[1, 1],
        macro_pars$Theta[2, 1]
      ),
      3,
      1
    )
  }
  
  dag_bar0
  
}

# Function to set initial microevolutionary parameter values
get_micro_init = function(dag_bar0, micro_pars, constraints, x, model) {
  list(
    dag_bar = dag_bar0,
    G = get_msb(
      micro_pars$phi_a,
      micro_pars$omega,
      micro_pars$M,
      constraints$beta,
      x = x,
      model = model
    )
  )
}

# Custom rOU function from POUMM so I don't have to install library on HPC
# rOU <- function(n, z0, t, alpha, theta, sigma) 
# {
#   ett <- exp(-alpha * t)
#   sd <- sigma * sqrt(t)
#   a <- alpha > 0
#   sd[a] <- sigma[a] * sqrt((1 - ett[a]^2)/(2 * alpha[a]))
#   mean <- z0 * ett + theta * (1 - ett)
#   nan <- is.infinite(t) & alpha == 0
#   mean[nan] <- z0[nan]
#   rnorm(n, mean = mean, sd = sd)
# }

# Shortcut function for calculating var_z
var_z <- function(var_ds, var_as, cov_da, B) {
  var_ds + B ^ 2 * var_as + 2 * B * cov_da
}

# Calculations of gsmax follow Sack and Buckley 2016
biophysical_constant <- function(D_wv, v) D_wv / v

morphological_constant <- function(c, h, j) {
  (pi * c ^ 2) / (j ^ 0.5 * (4 * h * j + pi * c))
}

# Faster mean and variance functions
mean1 <- function(x) sum(x) / length(x)
var1 <- function(x) sum((x - sum(x) / length(x)) ^ 2) / (length(x) - 1)

# Acceptable names in each parameter set
par_names <- function(type) {
  
  assert_character(type, len = 1L)
  type <- match.arg(type, c("ou", "qg", "trait", "physical_limits"))
  switch(
    type,
    ou = c("alpha", "sigma", "theta"),
    qg = c("N", "V_e", "V_mX", "V_mY", "r_XY", "V_s", "W_max"),
    trait = c("X_0", "Y_0", "z_0", "beta"),
    physical_limits = c("min_as", "max_as", "min_ds", "max_ds", "max_fs")
  )
  
}

# Check parameter inputs for simulate_lineage()
check_pars = function(ou_pars, qg_pars, trait_pars, physical_limits) {
  
  # Checks
  assert_list(ou_pars)
  assert_names(names(ou_pars), permutation.of = par_names("ou"))
  assert_number(ou_pars$alpha, lower = 0)
  assert_number(ou_pars$sigma, lower = 0)
  assert_number(ou_pars$theta)
  
  assert_list(qg_pars)
  assert_names(names(qg_pars), permutation.of = par_names("qg"))
  assert_number(qg_pars$N, lower = 0)
  assert_number(qg_pars$V_e, lower = 0)
  assert_number(qg_pars$V_mX, lower = 0)
  assert_number(qg_pars$V_mY, lower = 0)
  assert_number(qg_pars$r_XY, lower = -1, upper = 1)
  assert_number(qg_pars$V_s, lower = 0)
  assert_number(qg_pars$W_max, lower = 0)
  
  assert_list(trait_pars)
  assert_names(names(trait_pars), permutation.of = par_names("trait"))
  assert_number(trait_pars$X_0, null.ok = TRUE)
  assert_number(trait_pars$Y_0, null.ok = TRUE)
  assert_number(trait_pars$z_0)
  assert_number(trait_pars$beta)
  
  assert_list(physical_limits)
  assert_names(names(physical_limits), 
               permutation.of = par_names("physical_limits"))
  assert_number(physical_limits$min_as)
  assert_number(physical_limits$max_as)
  assert_number(physical_limits$min_ds)
  assert_number(physical_limits$max_ds)
  assert_number(physical_limits$max_fs)
  
}

# Calculate absolute fitness using univariate Gaussian
get_W <- function(z, W_max, u, V_s, log = FALSE) {
  if (log) {
    ret <- log(W_max) - (z - u) ^ 2 / V_s
  } else {
    ret <- W_max * exp(-(z - u) ^ 2 / V_s) 
  }
  ret
}

# Function to concatenate simulations from individual output files
concatenate_sims <- function(.x, output_name = NULL) {
  
  if (is.null(output_name)) output_name <- .x[1]
  
  output <- str_c("simulations/", .x)
  sims <- map(output, ~{
    str_c(.x, "/", list.files(.x))
  }) %>%
    unlist()
  
  pb <- progress_bar$new(length(sims),
                         format = "[:bar] :current/:total (:percent) :eta")
  
  sim_set <- seq_len(length(sims)) %>%
    map_dfr(~{
      pb$tick()
      read_rds(sims[.x])
    })
  
  write_rds(sim_set, glue("objects/{output_name}.rds"))
  
}

# Functions to summarize scaling estimates from simulations
est_scaling = function(.x, explanatory, response) {
  
  .f = glue("{response} ~ {explanatory}")
  ols = lm(.f, data = .x) %>%
    tidy(conf.int = TRUE, conf.level = 0.99) %>%
    filter(term == explanatory) %>%
    select(term, estimate, conf.low, conf.high) %>%
    mutate(across(where(is.numeric), ~ {-.x}), method = "ols")
  
  sma = sma(.f, data = .x, method = "SMA", alpha = 0.01) %>%
    extract2("coef") %>%
    extract2(1) %>%
    rownames_to_column() %>%
    mutate(term = "mean_Y", estimate = `coef(SMA)`, 
           conf.low = `lower limit`, conf.high = `upper limit`) %>%
    filter(rowname == "slope") %>%
    select(term, estimate, conf.low, conf.high) %>%
    mutate(across(where(is.numeric), ~ {-.x}), method = "sma")
  
  bind_rows(ols, sma)
  
}  

# Function to whether constituent traits satisfy physical limits
check_physical_limits = function(X, Y, physical_limits) {
  
  # Logical: does element exceed a physical limit?
  i = exp(X - log(1e6) + Y) > physical_limits$max_fs |
    X < physical_limits$min_ds |
    X > physical_limits$max_ds |
    Y < physical_limits$min_as |
    Y > physical_limits$max_as
  
  ret = any(i)
  attr(ret, "i") = which(i)

  ret
  
}

# Functions to calculate expected covariance between traits with truncated bivariate normal distribution


get_cov = function(Mu, Sigma, B, z_max) {
  
  assert_numeric(Mu, len = 2)
  assert_matrix(Sigma, nrows = 2, ncols = 2)
  assert_true(all(eigen(Sigma)$values > 0))
  assert_numeric(B, len = 2)
  assert_numeric(z_max)
  
  denominator = get_integral(Mu, Sigma, B, z_max, "denominator")
  
  E1_num = get_E_numerator(1, Mu, Sigma, B, z_max)
  E2_num = get_E_numerator(2, Mu, Sigma, B, z_max)
  
  E1 = E1_num / denominator
  E2 = E2_num / denominator

  E1_2_num = get_E_2_numerator(1, Mu, Sigma, B, z_max)
  E2_2_num = get_E_2_numerator(2, Mu, Sigma, B, z_max)
  
  E1_2 = E1_2_num / denominator
  E2_2 = E2_2_num / denominator
  
  # integrate to get cov
  c(
    mu1 = E1,
    mu2 = E2,
    var1 = E1_2 - E1 ^ 2,
    var2 = E2_2 - E2 ^ 2,
    cov = get_integral(Mu, Sigma, B, z_max, "C_numerator", E1, E2) / denominator 
  )
  
}

get_integral = function(Mu, Sigma, B, z_max, which, E1 = NULL, E2 = NULL) {
  which = switch(
    which,
    E_numerator = "E_numerator",
    E2_numerator = "E2_numerator",
    C_numerator = "C_numerator",
    denominator = "denominator"
  )
  assert_number(E1, null.ok = TRUE)
  assert_number(E2, null.ok = TRUE)
  integrate(
    f = function(x1, Mu, Sigma, B, z_max, which, E1, E2) {
      map_dbl(x1, \(x_fix) {
        integrate(
          f = function(x_var, x_fix, Mu, Sigma, B, z_max, which, E1, E2) {
            
            x1 = if (which == "E_numerator") {
              x_var
            } else if (which == "E2_numerator") {
              x_var ^ 2
            } else if (which == "C_numerator") {
              (x_var - E1) * (x_fix - E2)
            } else {
              1
            }
            
            x2 = matrix(c(x_var, rep(x_fix, length(x_var))), ncol = 2)
            
            x1 * exp(mvnfast::dmvn(x2, mu = Mu, sigma = Sigma, log = TRUE))
            
          },
          lower = -Inf,
          upper = (z_max - B[2] * x_fix) / B[1],
          x_fix = x_fix,
          Mu = Mu,
          Sigma = Sigma,
          B = B,
          z_max = z_max,
          which = which,
          E1 = E1,
          E2 = E2,
          subdivisions = 1000L,
          abs.tol = 1e-10
        )$value
      })
    },
    lower = -Inf,
    upper = Inf,
    Mu = Mu,
    Sigma = Sigma,
    B = B,
    z_max = z_max,
    which = which,
    E1 = E1,
    E2 = E2,
    subdivisions = 1000L,
    abs.tol = 1e-10
  )$value
}

get_E_numerator = function(n, Mu, Sigma, B, z_max) {
  
  assert_int(n, lower = 1, upper = 2)
  o = c(n, -n + 3)
  Sigma1 = Sigma
  diag(Sigma1) = diag(Sigma)[o]
  get_integral(Mu[o], Sigma1, B[o], z_max, "E_numerator")
  
}

get_E_2_numerator = function(n, Mu, Sigma, B, z_max) {
  
  assert_int(n, lower = 1, upper = 2)
  o = c(n, -n + 3)
  Sigma1 = Sigma
  diag(Sigma1) = diag(Sigma)[o]
  get_integral(Mu[o], Sigma1, B[o], z_max, "E2_numerator")
  
}

# Multivariate OU functions
# Functions to calculate covariance matrix of traits in multivariate OU process at time t
MVOU_integrand = function(s, Alpha, sigma, index) {
  map_dbl(s,
          \(.s, Alpha, sigma, index) {
            mat = expm(-Alpha * .s) %*% sigma %*% t(sigma) %*% expm(-t(Alpha) * .s)
            mat[index[1], index[2]]
          },
          Alpha = Alpha,
          sigma = sigma,
          index = index)
}

cov_MVOU1 = function(t, Alpha, Sigma, check = FALSE) {
  assert_flag(check)
  if (check) {
    assert_number(t, lower = 0)
    assert_matrix(Alpha)
    assert_matrix(Sigma)
    assert_true(nrow(Alpha) == ncol(Alpha))
    assert_true(all(dim(Alpha) == dim(Sigma)))
  }
  
  m = nrow(Alpha)
  
  crossing(i = seq_len(m), j = seq_len(m)) |>
    pmap_dbl(\(i, j, t, Alpha, sigma) {
      integrate(
        MVOU_integrand,
        lower = 0,
        upper = t,
        index = c(i, j),
        Alpha = Alpha,
        sigma = sigma
      )$value
    },
    t = t,
    Alpha = Alpha,
    sigma = sqrtm(Sigma)) |>
    matrix(nrow = m)
  
}

cov_MVOU = function(t, Alpha, Sigma, check = FALSE) {
  assert_flag(check)
  if (check) {
    assert_numeric(t, lower = 0)
    assert_matrix(Alpha)
    assert_matrix(Sigma)
    assert_true(nrow(Alpha) == ncol(Alpha))
    assert_true(all(dim(Alpha) == dim(Sigma)))
  }
  
  map(t, \(.t, Alpha, Sigma) {
    cov_MVOU1(.t, Alpha, Sigma, check = FALSE)
  }, Alpha = Alpha, Sigma = Sigma)
  
}

rMVOU = function(n, Z0, t, Alpha, Theta, Sigma) {
  
  k = nrow(Alpha)
  Mu = expm(-Alpha * t) %*% Z0 + (diag(k) - expm(-Alpha * t)) %*% Theta
  SigmaOU = cov_MVOU1(t, Alpha, Sigma)
  mvnfast::rmvn(n, Mu, SigmaOU)
  
}


# From mvMORPH:::StationaryVariance()
StationaryVariance = function(Alpha, Sigma)
{
  alpha = sqrtm(Alpha)
  Sigma <- Sigma
  eig <- eigen(alpha)
  P <- eig$vectors
  invP <- solve(P)
  eigvalues <- eig$values
  p = dim(Sigma)[1]
  Mat <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      Mat[i, j] <- 1 / (eigvalues[i] + eigvalues[j])
    }
  }
  StVar <- P %*% (Mat * (invP %*% Sigma %*% t(invP))) %*% t(P)
  return(StVar)
}

par2mat = function(par, r_alpha, r_sigma) {
  A_Y = exp(par[1])
  A_Z = exp(par[2])
  S_Y = exp(par[3])
  S_Z = exp(par[4])
  A = matrix(c(A_Y, rep(r_alpha * sqrt(A_Y * A_Z), 2), A_Z), 2, 2)
  S = matrix(c(S_Y, rep(r_sigma * sqrt(S_Y * S_Z), 2), S_Z), 2, 2)
  list(Alpha = A, Sigma = S)
}

f1 = function(par, r_alpha, r_sigma, S0) {
  mat = par2mat(par, r_alpha, r_sigma)
  S1 = StationaryVariance(mat$Alpha, mat$Sigma)
  sum((diag(S1) - diag(S0)) ^ 2)
}

# Function to modify Alpha and Sigma matrices to keep stationary variance of 
# Y and Z constant, but change the covariance between Y and Z
update_macro_pars = function(macro_pars, r_alpha, r_sigma) {
  
  macro_pars$Alpha0 = macro_pars$Alpha
  macro_pars$Sigma0 = macro_pars$Sigma
  S0 = StationaryVariance(macro_pars$Alpha0, macro_pars$Sigma0)
  
  # initial guess  
  par = log(c(
    macro_pars$Alpha0[1, 1],
    macro_pars$Alpha0[2, 2],
    macro_pars$Sigma0[1, 1],
    macro_pars$Sigma0[2, 2]
  ))
  
  fit = optim(par, f1, r_alpha = r_alpha, r_sigma = r_sigma, S0 = S0)
  mat = par2mat(fit$par, r_alpha, r_sigma)
  
  macro_pars$Alpha = mat$Alpha
  macro_pars$Sigma = mat$Sigma
  macro_pars
  
}

rTrajectoryMVOU = function(t, Z0, delta_t, Alpha, Theta, Sigma, seed) {
  if (is.null(seed)) seed = as.integer(sample.int(1e9, 1))
  # use_python("~/.pyenv/shims/python") # laptop - need to install 3.12.0 using pyenv on other machines
  # use_python("/opt/anaconda3/bin/python")
  source_python("python/rTrajectoryMVOU.py")
  np = import("numpy", convert = FALSE)
  np$array(rTrajectoryMVOU_py(t, Z0, delta_t, sqrtm(Alpha), Theta, Sigma, seed)) |>
    py_to_r()
}

# G* for hypotheses 2-4
get_H = function(..., model) {
  switch(
    model,
    "h1" = get_H_h1(...),
    "h2" = get_H_h2(...),
    "h3" = get_H_h3(...)
  )
}

get_H_h1 = function(omega, beta) {
  matrix(c(-1 / omega, rep(-beta / omega, 2), beta^2 / omega), 2, 2)
}

get_H_h2 = function(phi_f, omega, beta) {
  matrix(c(
    -(1 + phi_f) / (omega * phi_f),
    rep(-(1 + beta * phi_f) / (omega * phi_f), 2),
    -(1 + beta^2 * phi_f) / (omega * phi_f)
  ), 2, 2)
}

get_H_h3 = function(phi_a, omega, beta) {
  matrix(c(
    -1 / omega,
    rep(-beta / omega, 2),
    -(phi_a * beta^2 + 1) / (phi_a * omega)
  ), 2, 2)
}

# Function to define the system of equations to solve for G-matrix as mutation-selection balance
# Python script showed that there may be no analytical solution or it is extremely complicated
# I used chatGPT to help with code for the numerical solution
solve_G = function(unknowns, H, M) {
  
  # Assign unknowns to G
  G = matrix(unknowns, nrow = 2, byrow = TRUE)
  
  # Compute GHG 
  # H is curvature of microevolutionary adaptive landscape (Lande 1980)
  GHG = G %*% H %*% G
  
  # Return the difference between GHG and M
  return(as.vector(GHG + M))
  
}

# Mutation-selection balance for G of both d_s and a_s
get_msb = function(..., x, model) {
  switch(
    model,
    "h2" = get_msb_h2(..., x = x),
    "h3" = get_msb_h3(..., x = x)
  )
}

get_msb_h2 = function(phi_f, omega, M, beta, x = 1e1 * M) {
  
  assert_number(phi_f, lower = 0)
  assert_number(omega, lower = 0)
  assert_matrix(M, nrow = 2, ncol = 2)
  assert_number(beta)
  
  H = get_H(phi_f, omega, beta, model = "h2")
  
  solution = nleqslv(
    x = x,
    fn = solve_G,
    control = list(allowSingular = TRUE),
    H = H,
    M = M
  )
  
  # Fix sign of solution
  if (solution$x[1] < 0) {
    d1 = solve_G(solution$x, H = H, M = M)
    d2 = solve_G(-solution$x, H = H, M = M)
    if (all(d1 == d2)) {
      solution$x = -solution$x
    }
  }
  
  if (solution$message == "Function criterion near zero") {
    return(matrix(solution$x, 2, 2))
  } else {
    warning("No solution found")
    return(NA)
  }
  
  
}

get_msb_h3 = function(phi_a, omega, M, beta, x = 1e1 * M) {
  
  assert_number(phi_a, lower = 0)
  assert_number(omega, lower = 0)
  assert_matrix(M, nrow = 2, ncol = 2)
  assert_number(beta)
  
  H = get_H(phi_a, omega, beta, model = "h3")
  
  solution = nleqslv(
    x = x,
    fn = solve_G,
    control = list(allowSingular = TRUE, xtol = 1e-16),
    H = H,
    M = M
  )

  # Fix sign of solution
  if (solution$x[1] < 0) {
    d1 = solve_G(solution$x, H = H, M = M)
    d2 = solve_G(-solution$x, H = H, M = M)
    if (all(d1 == d2)) {
      solution$x = -solution$x
    }
  }

  if (solution$message == "Function criterion near zero") {
    return(matrix(solution$x, 2, 2))
  } else {
    warning("No solution found")
    return(NA)
  }

}

# Functions to derive mu*
get_mu_star = function(..., model) {
  
  # Initial guess
  init_pars = rep(0, 2)
  
  solution = nleqslv(
    x = init_pars,
    fn = solve_mu,
    control = list(allowSingular = TRUE),
    ...,
    model = model
  )
  
  solution$x
  
}

solve_mu = function(unknowns, ..., model) {
  
  # Assign unknowns to mu0 (interspecific means of d_s and a_s)
  mu0 = unknowns
  
  # Compute stationary distribution
  mu1 = c(
    get_mu_dS1(mu0 = mu0, ..., model = model),
    get_mu_aS1(mu0 = mu0, ..., model = model)
  )
  
  # Return the difference between XY and XY1
  return(as.vector(mu1 - unknowns))
  
}

get_mu_dS1 = function(mu0, ..., model) {
  
  switch(
    model,
    "h2" = mu_dS1_h2(mu0, ...),
    "h3" = mu_dS1_h3(mu0, ...)
  )
  
}

get_mu_aS1 = function(mu0, ..., model) {
  
  switch(
    model,
    "h2" = mu_aS1_h2(mu0, ...),
    "h3" = mu_aS1_h3(mu0, ...)
  )
  
}

# mu_dS after 1 generation of selection in h2
mu_dS1_h2 = function(mu0, G, omega, phi_f, z0_gmax, z0_fS, beta, fS_min, g_opt) {
  
  dS = mu0[1]
  aS = mu0[2]
  Gd = G[1,1]
  Gda = G[1,2]
  
  a0 = fS_min * (Gd + Gda) / (omega * phi_f) -
    (Gd * (z0_gmax + z0_fS / phi_f) / omega + Gda * (beta * z0_gmax + z0_fS / phi_f) / omega)
  a1 = (1 - Gd * (1 + 1/phi_f) / omega - Gda * (beta + 1 / phi_f) / omega)
  a2 = -(Gd * (beta + 1 / phi_f) / omega + Gda * (beta ^ 2 + 1 / phi_f) / omega)
  a3 = (Gd + beta * Gda) / omega
  
  a0 + a1 * dS + a2 * aS + a3 * g_opt
  
}

# mu_aS after 1 generation of selection in h2
mu_aS1_h2 = function(mu0, G, omega, phi_f, z0_gmax, z0_fS, beta, fS_min, g_opt) {
  
  dS = mu0[1]
  aS = mu0[2]
  Ga = G[2,2]
  Gda = G[1,2]
  
  a0 = - Ga * (beta * z0_gmax - fS_min / phi_f + z0_fS / phi_f) / omega - Gda * (z0_gmax - fS_min / phi_f + z0_fS / phi_f) / omega
  a1 = - Ga * (beta + 1 / phi_f) / omega - Gda * (1 + 1/phi_f)/ omega
  a2 = 1 - Ga * (beta^2 + 1 / phi_f) / omega - Gda * (beta + 1 / phi_f) / omega
  a3 = (beta * Ga + Gda) / omega
  
  a0 + a1 * dS + a2 * aS + a3 * g_opt
  
}


# Functions to derive V*
get_V_star = function(x, ..., model) {
  
  solution = nleqslv(
    x = x,
    fn = solve_V,
    control = list(allowSingular = TRUE, xtol = 1e-16),
    ...,
    model = model
  )
  
  matrix(rep(solution$x, c(1, 2, 1)), 2, 2)
  
}

solve_V = function(unknowns, ..., model) {
  
  # Assign unknowns to V0 (interspecific covariance of d_s and a_s)
  V0 = matrix(c(unknowns[1], rep(unknowns[2], 2), unknowns[3]), nrow = 2, byrow = TRUE)
  
  # Compute stationary distribution
  V1 = c(
    get_Vd(V0 = V0, ..., model = model),
    get_Vda(V0 = V0, ..., model = model),
    get_Va(V0 = V0, ..., model = model)
  )
  
  # Return the difference between XY and XY1
  return(as.vector(abs(V1 - unknowns)))
  
}

# Variance in d_s
get_Vd = function(V0, ..., model) {
  switch(
    model,
    "h1" = get_Vd_h1(V0, ...),
    "h2" = get_Vd_h2(V0, ...),
    "h3" = get_Vd_h3(V0, ...)
  )
}

# Covariance of d_s and a_s
get_Vda = function(V0, ..., model) {
  switch(
    model,
    "h1" = get_Vda_h1(V0, ...),
    "h2" = get_Vda_h2(V0, ...),
    "h3" = get_Vda_h3(V0, ...)
  )
}

# Variance in a_s
get_Va = function(V0, ..., model) {
  switch(
    model,
    "h1" = get_Va_h1(V0, ...),
    "h2" = get_Va_h2(V0, ...),
    "h3" = get_Va_h3(V0, ...)
  )
}

# H1 functions
get_Vd_h1 = function(V0, V_gopt, G, beta, omega) {
  a1 = 1 - (beta * G[1, 2] + G[1, 1]) / omega
  a2 = (- beta ^ 2 * G[1, 2] - beta * G[1, 1]) / omega
  a3 = (beta * G[1, 2] + G[1, 1]) / omega
  # a1 ^ 2 * V0[1,1] + a2 ^ 2 * V0[2,2] + a3 ^ 2 * V_gopt + 2 * a1 * a2 * V0[1,2]
  V_dS_gopt = V_gopt / 2 + V0[1,1] / 2 - beta * V0[2,2] / 2
  V_aS_gopt = (V_gopt + beta ^ 2 * V0[2,2] - V0[1,1]) / (2 * beta)
  
  a1 ^ 2 * V0[1,1] + a2 ^ 2 * V0[2,2] + a3 ^ 2 * V_gopt + 
    2 * a1 * a2 * V0[1,2] + 2 * a1 * a3 * V_dS_gopt + 2 * a2 * a3 * V_aS_gopt
  
}

get_Va_h1 = function(V0, V_gopt, G, beta, omega) {
  a1 = (-beta * G[2, 2] - G[1, 2]) / omega
  a2 = (1 - (beta ^ 2 * G[2, 2] + G[1, 2] * beta) / omega)
  a3 = (beta * G[2, 2] + G[1, 2]) / omega
  # a1 ^ 2 * V0[1,1] + a2 ^ 2 * V0[2,2] + a3 ^ 2 * V_gopt + 2 * a1 * a2 * V0[1,2]
  V_dS_gopt = V_gopt / 2 + V0[1,1] / 2 - beta * V0[2,2] / 2
  V_aS_gopt = (V_gopt + beta ^ 2 * V0[2,2] - V0[1,1]) / (2 * beta)
  
  a1 ^ 2 * V0[1,1] + a2 ^ 2 * V0[2,2] + a3 ^ 2 * V_gopt + 
    2 * a1 * a2 * V0[1,2] + 2 * a1 * a3 * V_dS_gopt + 2 * a2 * a3 * V_aS_gopt
}

get_Vda_h1 = function(V0, V_gopt, G, beta, omega) {
  
  # Derivation is based on bilinearity property of covariance (with help from ChatGPT)
  # Cov(aX+bY,cW+dZ)=ac⋅Cov(X,W)+ad⋅Cov(X,Z)+bc⋅Cov(Y,W)+bd⋅Cov(Y,Z)
  
  a1x = 1 - (beta * G[1, 2] + G[1, 1]) / omega
  a2x = (- beta ^ 2 * G[1, 2] - beta * G[1, 1]) / omega
  a3x = (beta * G[1, 2] + G[1, 1]) / omega
  
  a1y = (- beta * G[2, 2] - G[1, 2]) / omega
  a2y = (1 - (beta ^ 2 * G[2, 2] + G[1, 2] * beta) / omega)
  a3y = (beta * G[2, 2] + G[1, 2]) / omega
  
  # a1x * a1y * V0[1,1] + (a1x * a2y + a2x * a1y) * V0[1,2] + a2x * a2y * V0[2,2] + a3x * a3y * V_gopt
  V_dS_gopt = V_gopt / 2 + V0[1,1] / 2 - beta * V0[2,2] / 2
  V_aS_gopt = (V_gopt + beta ^ 2 * V0[2,2] - V0[1,1]) / (2 * beta)
  
  a1x * a1y * V0[1,1] + 
    a2x * a2y * V0[2,2] + 
    a3x * a3y * V_gopt +
    (a1x * a2y + a2x * a1y) * V0[1,2] + 
    (a1x * a3y + a3x * a1y) * V_dS_gopt +
    (a2x * a3y + a3x * a2y) * V_aS_gopt
  
}

# H2 functions
get_Vd_h2 = function(V0, V_gopt, G, beta, omega, phi_f) {
  a1 = (1 - G[1,1] * (1 + 1/phi_f) / omega - G[1,2] * (beta + 1 / phi_f) / omega)
  a2 = -(G[1,1] * (beta + 1 / phi_f) / omega + G[1,2] * (beta ^ 2 + 1 / phi_f) / omega)
  a3 = (G[1,1] + beta * G[1,2]) / omega

  V_dS_gopt = V_gopt / 2 + V0[1,1] / 2 - beta * V0[2,2] / 2
  V_aS_gopt = (V_gopt + beta ^ 2 * V0[2,2] - V0[1,1]) / (2 * beta)

  a1 ^ 2 * V0[1,1] + a2 ^ 2 * V0[2,2] + a3 ^ 2 * V_gopt + 
    2 * a1 * a2 * V0[1,2] + 2 * a1 * a3 * V_dS_gopt + 2 * a2 * a3 * V_aS_gopt
  
}

get_Va_h2 = function(V0, V_gopt, G, beta, omega, phi_f) {
  a1 = - G[2,2] * (beta + 1 / phi_f) / omega - G[1,2] * (1 + 1/phi_f)/ omega
  a2 = 1 - G[2,2] * (beta^2 + 1 / phi_f) / omega - G[1,2] * (beta + 1 / phi_f) / omega
  a3 = (beta * G[2,2] + G[1,2]) / omega
  
  V_dS_gopt = V_gopt / 2 + V0[1,1] / 2 - beta * V0[2,2] / 2
  V_aS_gopt = (V_gopt + beta ^ 2 * V0[2,2] - V0[1,1]) / (2 * beta)
  
  a1 ^ 2 * V0[1,1] + a2 ^ 2 * V0[2,2] + a3 ^ 2 * V_gopt + 
    2 * a1 * a2 * V0[1,2] + 2 * a1 * a3 * V_dS_gopt + 2 * a2 * a3 * V_aS_gopt
    
}

get_Vda_h2 = function(V0, V_gopt, G, beta, omega, phi_f) {
  
  a1x = (1 - G[1,1] * (1 + 1/phi_f) / omega - G[1,2] * (beta + 1 / phi_f) / omega)
  a2x = -(G[1,1] * (beta + 1 / phi_f) / omega + G[1,2] * (beta ^ 2 + 1 / phi_f) / omega)
  a3x = (G[1,1] + beta * G[1,2]) / omega
  
  a1y = - G[2,2] * (beta + 1 / phi_f) / omega - G[1,2] * (1 + 1/phi_f)/ omega
  a2y = 1 - G[2,2] * (beta^2 + 1 / phi_f) / omega - G[1,2] * (beta + 1 / phi_f) / omega
  a3y = (beta * G[2,2] + G[1,2]) / omega
  
  # Derivation and test. can delete after chekcing
  # Cov(dS, g_opt)
  # g_opt = dS + B * aS
  # aS = (g_opt - dS) / B
  # Var(aS) = (1/B)^2 * Var_gopt + (1/B)^2 * Var_dS - 2 / B ^ 2 * Cov(gopt, dS)
  # Cov(gopt, dS) = ((1/B)^2 * Var_gopt + (1/B)^2 * Var_dS - Var(aS)) * B ^ 2 / 2
  # Cov(gopt, dS) = Var_gopt / 2 + Var_dS / 2 - B^ 2 * Var(aS) / 2

  # Cov(aS, g_opt)
  # g_opt = dS + B * aS
  # dS = g_opt - B * aS
  # Var(dS) = var(g_opt) + B^2 * var(aS) - 2 * B * cov(g_opt, aS)
  # Cov(g_opt, aS) = (var(g_opt) + B^2 * var(aS) - Var(dS)) / (2 * B)

  # dS = rnorm(100)
  # aS = rnorm(100)
  # B = 0.4
  # g_opt = dS + B * aS
  # var(aS)
  # (1/B)^2 * var(g_opt) + (1/B)^2 * var(dS) - 2 * 1 / B ^ 2 * cov(g_opt, dS)
  # var(g_opt)/ 2 + var(dS)/2 - B^ 2 * var(aS) / 2
  # (var(g_opt) + B^2 * var(aS) - var(dS)) / (2 * B)
  # cov(g_opt, aS)
  
  V_dS_gopt = V_gopt / 2 + V0[1,1] / 2 - beta * V0[2,2] / 2
  V_aS_gopt = (V_gopt + beta ^ 2 * V0[2,2] - V0[1,1]) / (2 * beta)
  
  a1x * a1y * V0[1,1] + 
    (a1x * a2y + a2x * a1y) * V0[1,2] + 
    (a1x * a3y + a3x * a1y) * V_dS_gopt +
    a2x * a2y * V0[2,2] + 
    (a2x * a3y + a3x * a2y) * V_aS_gopt +
    a3x * a3y * V_gopt
  
}

# Functions to summarize simulation variances and bootstrap CIs
boot_var = function(x, n_boot = 1e3, probs = c(0.025, 0.975)) {
  var_name = as_string(ensym(x))
  var_name_estimate = str_c("V_", var_name, "_estimate")
  var_name_lower = str_c("V_", var_name, "_lower")
  var_name_upper = str_c("V_", var_name, "_upper")
  q = quantile(replicate(n_boot, var(sample(x, length(
    x
  ), TRUE))), probs = probs)
  tibble(!!var_name_estimate := var(x),
         !!var_name_lower := q[1],
         !!var_name_upper := q[2])
}

boot_cov = function(x, y, n_boot = 1e3, probs = c(0.025, 0.975)) {
  var_x_name = as_string(ensym(x))
  var_y_name = as_string(ensym(y))
  var_name = str_c("V_", var_x_name, var_y_name)
  var_name_estimate = str_c(var_name, "_estimate")
  var_name_lower = str_c(var_name, "_lower")
  var_name_upper = str_c(var_name, "_upper")
  q = quantile(replicate(n_boot, {
    i = sample.int(length(x), length(x), TRUE)
    cov(x[i], y[i])
  }), probs = probs)
  tibble(!!var_name_estimate := cov(x,y),
         !!var_name_lower := q[1],
         !!var_name_upper := q[2])
}

# Wrapper to try different starting values when getting G*
redo = function(G_star) {
  if (is.na(G_star[1])) {
    return(TRUE)
  } else {
    if (G_star[1, 1] < 0 | G_star[2, 2] < 0) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
}

get_Gstar = function(phi, omega, M, beta, model) {
  G_star = get_msb(phi, omega, M, beta, x = 1e2 * M, model = model)
  if (redo(G_star)) {
    G_star = get_msb(phi, omega, M, beta, x = 1e1 * M, model = model)
  }
  if (redo(G_star)) {
    G_star = get_msb(phi, omega, M, beta, x = 1e3 * M, model = model)
  }
  G_star
}
