rm(list = ls())

source("r/functions.R")

library(ape)
library(broom)
library(checkmate)
library(cowplot)
library(ellipse)
library(expm)
library(furrr)
library(future)
library(glue)
library(magrittr)
library(mvnfast)
library(mvSLOUCH)
library(nleqslv)
library(PCMBaseCpp)
library(phylolm)
library(phytools)
library(POUMM)
library(progress)
library(RefManageR)
library(reticulate)
library(rlang)
library(scales)
library(taxonlookup)
library(tidybayes)
library(tidyverse)
library(tikzDevice)
library(units)

theme_set(theme_cowplot())

n_workers = parallel::detectCores(logical = FALSE)
n_spp = 1e3 # number of species for simulations

hypotheses = c(
  "Stomatal-area adaptation",
  "Stomatal-area minimization",
  "Stomatal adaptation + bounded size"
)
