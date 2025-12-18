# Package index

## Bayesian Statistics

Functions for Bayesian statistics.

- [`conjugate()`](https://danforthcenter.github.io/pcvr/reference/conjugate.md)
  : Bayesian testing using conjugate priors and method of moments for
  single or multi value traits.
- [`growthSS()`](https://danforthcenter.github.io/pcvr/reference/growthSS.md)
  : Ease of use growth model helper function.
- [`fitGrowth()`](https://danforthcenter.github.io/pcvr/reference/fitGrowth.md)
  : Ease of use wrapper function for fitting various growth models
  specified by growthSS
- [`growthPlot()`](https://danforthcenter.github.io/pcvr/reference/growthPlot.md)
  : Function to visualize models made by fitGrowth.
- [`stat_brms_model()`](https://danforthcenter.github.io/pcvr/reference/stat_growthss.md)
  [`stat_growthss()`](https://danforthcenter.github.io/pcvr/reference/stat_growthss.md)
  [`stat_nlme_model()`](https://danforthcenter.github.io/pcvr/reference/stat_growthss.md)
  [`stat_lme_model()`](https://danforthcenter.github.io/pcvr/reference/stat_growthss.md)
  [`stat_nlrq_model()`](https://danforthcenter.github.io/pcvr/reference/stat_growthss.md)
  [`stat_nls_model()`](https://danforthcenter.github.io/pcvr/reference/stat_growthss.md)
  : Show brms model in ggplot layer
- [`testGrowth()`](https://danforthcenter.github.io/pcvr/reference/testGrowth.md)
  : Hypothesis testing for fitGrowth models.
- [`barg()`](https://danforthcenter.github.io/pcvr/reference/barg.md) :
  Function to help fulfill elements of the Bayesian Analysis Reporting
  Guidelines.
- [`combineDraws()`](https://danforthcenter.github.io/pcvr/reference/combineDraws.md)
  : Combine Draws From brms Models
- [`plotPrior()`](https://danforthcenter.github.io/pcvr/reference/plotPrior.md)
  : Check priors used in ease of use brms functions
- [`distributionPlot()`](https://danforthcenter.github.io/pcvr/reference/distributionPlot.md)
  : Function for plotting iterations of posterior distributions
- [`brmViolin()`](https://danforthcenter.github.io/pcvr/reference/brmViolin.md)
  : Function to visualize hypotheses tested on brms models similar to
  those made using growthSS outputs.

## Modeling

Functions for Bayesian or Frequentist Longitudinal (and other) Modeling

- [`growthSim()`](https://danforthcenter.github.io/pcvr/reference/growthSim.md)
  : Growth data simulating function
- [`growthSS()`](https://danforthcenter.github.io/pcvr/reference/growthSS.md)
  : Ease of use growth model helper function.
- [`fitGrowth()`](https://danforthcenter.github.io/pcvr/reference/fitGrowth.md)
  : Ease of use wrapper function for fitting various growth models
  specified by growthSS
- [`growthPlot()`](https://danforthcenter.github.io/pcvr/reference/growthPlot.md)
  : Function to visualize models made by fitGrowth.
- [`stat_brms_model()`](https://danforthcenter.github.io/pcvr/reference/stat_growthss.md)
  [`stat_growthss()`](https://danforthcenter.github.io/pcvr/reference/stat_growthss.md)
  [`stat_nlme_model()`](https://danforthcenter.github.io/pcvr/reference/stat_growthss.md)
  [`stat_lme_model()`](https://danforthcenter.github.io/pcvr/reference/stat_growthss.md)
  [`stat_nlrq_model()`](https://danforthcenter.github.io/pcvr/reference/stat_growthss.md)
  [`stat_nls_model()`](https://danforthcenter.github.io/pcvr/reference/stat_growthss.md)
  : Show brms model in ggplot layer
- [`testGrowth()`](https://danforthcenter.github.io/pcvr/reference/testGrowth.md)
  : Hypothesis testing for fitGrowth models.

## Lemnatech System Support

Functions for use with Lemnatech style longitudinal data.

- [`pcv.outliers()`](https://danforthcenter.github.io/pcvr/reference/pcv.outliers.md)
  : Remove outliers from bellwether data using cook's distance
- [`pcv.time()`](https://danforthcenter.github.io/pcvr/reference/pcv.time.md)
  : Time conversion and plotting for bellwether data
- [`pcv.water()`](https://danforthcenter.github.io/pcvr/reference/pcv.water.md)
  : Read in lemnatech watering data from metadata.json files
- [`pwue()`](https://danforthcenter.github.io/pcvr/reference/pwue.md) :
  Calculate pseudo water use efficiency from phenotype and watering data

## Single-Value Traits

Functions suited for use with single-value trait data (Height, Area,
etc)

- [`conjugate()`](https://danforthcenter.github.io/pcvr/reference/conjugate.md)
  : Bayesian testing using conjugate priors and method of moments for
  single or multi value traits.
- [`growthSS()`](https://danforthcenter.github.io/pcvr/reference/growthSS.md)
  : Ease of use growth model helper function.
- [`frem()`](https://danforthcenter.github.io/pcvr/reference/frem.md) :
  Variance partitioning using Full Random Effects Models
- [`relativeTolerance()`](https://danforthcenter.github.io/pcvr/reference/relativeTolerance.md)
  : Calculate relative tolerance of some phenotype(s) relative to
  control
- [`cumulativePheno()`](https://danforthcenter.github.io/pcvr/reference/cumulativePheno.md)
  : Reduce phenotypes in longitudinal data to cumulative sums of
  phenotypes.

## Multi-Value Traits

Functions suited for use with multi-value trait data (spectral indices,
colorspaces, CropReporter, etc)

- [`conjugate()`](https://danforthcenter.github.io/pcvr/reference/conjugate.md)
  : Bayesian testing using conjugate priors and method of moments for
  single or multi value traits.
- [`mv_ag()`](https://danforthcenter.github.io/pcvr/reference/mv_ag.md)
  : Multi Value Trait Aggregation function
- [`pcv.emd()`](https://danforthcenter.github.io/pcvr/reference/pcv.emd.md)
  [`pcv.euc()`](https://danforthcenter.github.io/pcvr/reference/pcv.emd.md)
  : Earth Mover's Distance between spectral histograms
- [`pcv.joyplot()`](https://danforthcenter.github.io/pcvr/reference/pcv.joyplot.md)
  : Make Joyplots for multi value trait plantCV data
- [`pcadf()`](https://danforthcenter.github.io/pcvr/reference/pcadf.md)
  : Function to run a PCA, plot and optionally return the data with PCA
  coordinates and pca object
- [`mvSS()`](https://danforthcenter.github.io/pcvr/reference/mvSS.md) :
  Ease of use multi-value trait model helper function.
- [`pcv.net()`](https://danforthcenter.github.io/pcvr/reference/pcv.net.md)
  : Network analysis of a distance matrix
- [`pcv.plsr()`](https://danforthcenter.github.io/pcvr/reference/pcv.plsr.md)
  : Run Partial Least Squares Regression on spectral data

## Reading PlantCV output

Functions for reading different formats of PlantCV output (note for
PlantCV \>= 4.1) data.table::fread or utils::read.csv are perfectly
adequate.

- [`read.pcv.3()`](https://danforthcenter.github.io/pcvr/reference/read.pcv.3.md)
  : Read in plantCV csv from bellwether phenotyper style experiments
  analyzed with plantCV versions \<4.
- [`read.pcv()`](https://danforthcenter.github.io/pcvr/reference/read.pcv.md)
  : Read in plantCV csv output in wide or long format
