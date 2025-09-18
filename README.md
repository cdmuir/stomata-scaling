[![DOI](https://zenodo.org/badge/1050530864.svg)](https://doi.org/10.5281/zenodo.17055962)

# stomata-scaling

The repository contains source code for analyses and simulations related to:

*Bounds on stomatal size can explain scaling with stomatal density in forest plants*

Congcong Liu, Christopher D. Muir, Lawren Sack, Ying Li, Jiahui Zhang, Guirui Yu, Li Xu, Mingxu Li, Zihao Zhang, Hugo Jan de Boer, Xingguo Han, Nianpeng He. Accepted in *New Phytologist*.

## Downloading repository 

1. Download or clone this repository to your machine.

```
git clone git@github.com:cdmuir/stomatal-scaling.git
```

2. Open `stomata-scaling.Rproj` in [RStudio](https://www.rstudio.com/)
3. Install *R* packages if necessary using `r/install-packages.R`

## Brief description of *R* scripts in `r/` directory

- `functions.R`: custom functions
- `header.R`: script to clear workspace and load packages

- `01_fit-phylolm.r`: fit **phylolm** model to estimate scaling exponent
- `02_plot-hypotheses.R`: make figure 1
- `03_plot-results.R`: make figure 3 showing results from `01_fit-phylolm.r`
- `04_fit-OUCH.R`: estimate OU covariance between stomatal size and density
- `05_sim-OUCH.R`: simulate synthetic data from OUCH model for parametric bootstrap
- `06_bootstrap-OUCH.R`: fit synthetic data from OUCH model for parameteric bootstrap
- `07_summarize-OUCH.R`: summarize OUCH parameter bootstrap (Table 3)
- `08_plot-ellipses.R`: make figure 4
- `09_set-parameters.R`: define baseline parameter sets for simulations
- `10_get-h1V*.R`: derive among-species size-density covariance (hypothesis 1)
- `11_get-h2G*.R`: derive within-species genetic covariance matrix, figures S2-S4 (hypothesis 2)
- `12_get-h2mu*.R`: derive among-species trait means (hypothesis 2)
- `13_get-h2V*.R`: derive among-species size-density covariance (hypothesis 2)
- `14_get-h3G*.R`: derive within-species genetic covariance matrix, figures S5-S7 (hypothesis 3)
- `15_compare-vectors.R`: make selection response vectors for figures S8-10
- `16_sim-h2-baseline.R`: Recursion equations to determine stationary variance of V* under baseline scenario (hypothesis 2)
- `17_sim-h3-baseline.R`: Recursion equations to determine stationary variance of V* under baseline scenario (hypothesis 3)
- `18_sim-h2-Mratio.R`: simulate effect of **M** ratio on V* (hypothesis 2)
- `19_sim-h3-Mratio.R`: simulate effect of **M** ratio on V* (hypothesis 3)
- `20_sim-h2-rM.R`: simulate effect of mutational correlation on V* (hypothesis 2)
- `21_sim-h3-rM.R`: simulate effect of mutational correlation on V* (hypothesis 3)
- `22_sim-h2-omega.R`: simulate effect of $\omega$ (hypothesis 2)
- `23_sim-h3-omega.R`: simulate effect of $\omega$ (hypothesis 3)
- `24_sim-h2-phi.R`: simulate effect of $\phi_f$ on V* (hypothesis 2)
- `25_sim-h3-phi.R`: simulate effect of $\phi_a$ on V* (hypothesis 3)
- `26_summarize-sims.R`: summarize simulation results
- `27_plot-sims.R`: make figures S11-S15
- `28_analyze-covariance.R`: analysis for section S4

## Raw data

- `raw-data/stomatal-data.csv`: stomatal size and density data
  + species
  + genus
  + family
  + stomatal density [1 / mm^2]
  + Stomatal length/size [um]
  + plant functional group/1-tree/2-shrub/3-herb
- `raw-data/phylogeny.tre`: phylogeny of species in stomatal data (newick format)
