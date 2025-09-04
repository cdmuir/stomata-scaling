source("r/header.R")

phy = read.tree("raw-data/Phylogeny.tre")

dat = read_csv("raw-data/stomatal-data.csv")

# Fix error
dat$family[which(dat$species == "Selaginella_tamariscina")] = "Selaginellaceae"

species_list <- intersect(phy$tip.label, dat$species)

phy %<>% keep.tip(species_list)
dat %<>% 
  filter(species %in% species_list) %>%
  select(
    species,
    genus,
    family,
    D_S = `stomatal density`,
    L = `Stomatal length/size`,
    SPI = `Stomatal pore area index`,
    pfg = `plant functional group/1-tree/2-shrub/3-herb`
  ) %>%
  mutate(
    grass = family == "Poaceae",
    c = 0.5,
    h = 0.5,
    j = ifelse(grass, 0.125, 0.5), # no effect if j = 0.36 as in de Boer et al. 2016
    m = morphological_constant(c, j, h),
    b = biophysical_constant(2.49e-5, 2.24e-2),
    A_S = L ^ 2 * j,
    g_max = b * m * D_S * sqrt(A_S),
    f_S = drop_units(set_units(set_units(D_S, mm ^ -2), um ^ -2) * 
                     set_units(A_S, um ^ 2))
  ) %>%
  filter(
    !is.na(D_S),
    !is.na(A_S)
  )

dat_spp <- dat %>%
  mutate(
    d_s = log(D_S),
    log_gmax = log(g_max),
    log_spi = log(SPI),
    log_fs = log(f_S),
    a_s = log(A_S),
    shrub = pfg == 2,
    tree = pfg == 1
  )

dat %<>% 
  group_by(species) %>%
  summarize(
    d_s = mean(log(D_S), na.rm = TRUE),
    log_gmax = mean(log(g_max), na.rm = TRUE),
    log_spi = mean(log(SPI), na.rm = TRUE),
    log_fs = mean(log(f_S), na.rm = TRUE),
    a_s = mean(log(A_S), na.rm = TRUE),
    n = length(species),
    j = first(j),
    m = first(m),
    b = first(b),
    pfg = floor(median(pfg)),
    family = first(family),
    grass = first(grass),
    .groups = "drop"
  ) %>%
  mutate(
    shrub = pfg == 2,
    tree = pfg == 1
  ) 

phy %<>% keep.tip(dat$species)

write_rds(phy, "objects/phy.rds")

tax_all = plant_lookup() %>%
  select(-genus) %>%
  distinct()

dat %<>% left_join(tax_all, by = "family")  
dat_spp %<>% left_join(tax_all, by = "family")  

dat$group[dat$family == "Leguminosae"] <- "Angiosperms" 
dat$group[dat$family == "Asterales"] <- "Angiosperms" 
dat$group[dat$family == "Cornales"] <- "Angiosperms" 
dat$group[dat$family == "Compositae"] <- "Angiosperms" 
dat$group[dat$family == "Xanthorrhoeaceae"] <- "Angiosperms" 
dat$group[dat$family == "Polypodiales"] <- "Pteridophytes" 

stopifnot(all(!is.na(dat$group)))

dat_spp$group[dat_spp$family == "Leguminosae"] <- "Angiosperms" 
dat_spp$group[dat_spp$family == "Asterales"] <- "Angiosperms" 
dat_spp$group[dat_spp$family == "Cornales"] <- "Angiosperms" 
dat_spp$group[dat_spp$family == "Compositae"] <- "Angiosperms" 
dat_spp$group[dat_spp$family == "Xanthorrhoeaceae"] <- "Angiosperms" 
dat_spp$group[dat_spp$family == "Polypodiales"] <- "Pteridophytes" 

stopifnot(all(!is.na(dat_spp$group)))

write_rds(dat, "objects/dat.rds")
write_rds(dat_spp, "objects/dat_spp.rds")
  
# Fit phylolm ----

plan(multisession, workers = 10)
set.seed(96969690)
fit_phylolm1 <- phylolm(
  formula = d_s ~ a_s * group + grass + shrub + tree,
  data = column_to_rownames(dat, "species"), 
  phy = phy,
  model = "OUrandomRoot",
  upper.bound = 10, 
  boot = 1e3
)

write_rds(fit_phylolm1, "objects/fit_phylolm1.rds")

plan(multisession, workers = 10)
set.seed(743125984)
outliers <- c("Torreya_fargesii")
fit_phylolm2 <- phylolm(
  formula = d_s ~ a_s * group + grass + shrub + tree,
  data = column_to_rownames(filter(dat, !(species %in% outliers)), "species"), 
  phy = drop.tip(phy, outliers),
  model = "OUrandomRoot",
  upper.bound = 10,
  boot = 1e3
)

write_rds(fit_phylolm2, "objects/fit_phylolm2.rds")

# Fit lm for comparison ----
fit_lm1 <- lm(
  formula = d_s ~ a_s * group + grass + shrub + tree,
  data = column_to_rownames(dat, "species")
)

write_rds(fit_lm1, "objects/fit_lm1.rds")

fit_lm2 <- lm(
  formula = d_s ~ a_s * group + grass + shrub + tree,
  data = column_to_rownames(filter(dat, !(species %in% outliers)), "species"), 
)

write_rds(fit_lm2, "objects/fit_lm2.rds")
