# Estimate OU covariance between stomatal size and density
source("r/header.R")

dat = read_rds("objects/dat.rds")
phy = read_rds("objects/phy.rds")
fit_simmap = make.simmap(phy, setNames(dat$group, dat$species))

# Fit all data ----
mdat = dat |>
  select(d_s, a_s) |>
  as.matrix()
rownames(mdat) = dat$species

rgs = fit_simmap$maps |>
  map_chr(\(.x) names(.x)[which.max(.x)])

fit_ouch = ouchModel(phy, mdat, regimes = rgs, Atype = "Symmetric")

write_rds(fit_simmap, "objects/fit_simmap.rds")
write_rds(fit_ouch, "objects/fit_ouch.rds")
