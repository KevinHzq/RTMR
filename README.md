# RTMR

<!-- badges: start -->
[![R-CMD-check](https://github.com/KevinHzq/RTMR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/KevinHzq/RTMR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

RTMR is an open-source R implementation of Risk Terrain Modeling (RTM).
It replicates the statistical procedure of the commercial RTMDx Utility
(Heffner, 2013) and extends it for epidemiological applications with rate
(offset) models, adjustment covariates, and holdout validation.

RTM diagnoses which features of the environment (bars, transit stops,
vacant lots, ...) are associated with where events such as crimes
concentrate, estimates *how far* each feature's influence reaches, and
combines the selected factors into a map of relative risk.

## How it works

1. A grid of cells is laid over the study area; outcome events are counted
   per cell.
2. Each risk factor is *operationalized* into candidate binary variables at
   a ladder of spatial influences: proximity ("within *d* of a feature")
   and/or density ("inside a high-concentration cluster at bandwidth *d*",
   Epanechnikov kernel, binarized at mean + 2 SD).
3. A cross-validated, sign-constrained penalized Poisson regression culls
   the candidates, protecting against spurious correlations from testing
   many variables.
4. Bidirectional stepwise regression by BIC selects the final factors and
   their optimal spatial influences, under both Poisson and negative
   binomial distributions; the best BIC wins. At most one variable — one
   operationalization at one distance — may represent each factor.
5. Each selected factor's exponentiated coefficient is its relative risk
   value (RRV; an adjusted incidence rate ratio), and the fitted surface,
   rescaled to its minimum, gives each cell's relative risk score.

## Installation

```r
# install.packages("remotes")
remotes::install_github("KevinHzq/RTMR")
```

## Usage

```r
library(RTMR)

fit <- rtm(
  outcome = crimes,                       # sf point layer of events
  factors = list(
    bars  = bars,                         # sf point layers of risk factors
    parks = list(data = parks, operation = "proximity"), # polygons: proximity only
    ops   = list(data = ops, type = "protective")
  ),
  boundary = city,                        # study area polygon (projected CRS)
  cell_size = 250,                        # grid resolution (~half a block)
  block_length = 500                      # average block length
)

fit                  # selected factors, spatial influences, RRVs
plot(fit)            # risk terrain map
write_rtm(fit, "risk_map.tif")            # GeoTiff for GIS (needs terra)

validate_rtm(fit, next_period_events)     # holdout capture/PAI, ROC/AUC
```

Rate models (e.g. spatial case fatality: fatal overdoses per overdose
event) are supported through the `offset` argument, and pre-computed
cell-level covariates — candidate risk factors or always-included
adjustment terms (continuous, or categorical with a chosen reference
level) — through `covariates`. See the vignette for the full workflow and
interpretation guidance:

```r
vignette("RTMR")
```

## References

- Caplan, J. M., Kennedy, L. W., & Piza, E. L. (2013). *Risk Terrain
  Modeling Diagnostics Utility User Manual (Version 1.0)*. Newark, NJ:
  Rutgers Center on Public Security.
- Heffner, J. (2013). Statistics of the RTMDx Utility. In the above manual.
- Chainey, S., Tompson, L., & Uhlig, S. (2008). The utility of hotspot
  mapping for predicting spatial patterns of crime. *Security Journal*,
  21, 4–28.

## License

MIT © Kevin Hu
