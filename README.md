
<!-- README.md is generated from README.Rmd. Please edit that file -->

# joineRmeta

<!-- badges: start -->

[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/mesudell/joineRmeta?branch=master&svg=true)](https://ci.appveyor.com/project/mesudell/joineRmeta)
<!-- badges: end -->

The `joineRmeta` package implements methods to analyse multi-study joint
data consisting of a single continuous longitudinal outcome, and a
single possible censored time-to-event outcome. The modelling framework
for the longitudinal data is a linear mixed effects model (Laird and
Ware, 1982). The modelling framework for the time-to-event outcome is a
Cox proportional hazards model with an unspecified baseline hazard (Cox,
1972). The longitudinal and time-to-event sub-model are linked through
an association structure. Currently only the random effects only
proportional association is available (see Gould et al, 2015). The
methodology used to fit the model is described in Henderson et al (2000)
and Wulfsohn and Tsiatis (1997).

The `joineRmeta` package contains methods to perform the second stage of
a two stage meta-analysis (MA) of study specific joint modelling fits,
and one stage MA of multi-study joint data where between study
heterogeneity can be accounted for using interaction terms with study
membership variables, study level random effects, or baseline hazard
stratified by study. The package also contains plotting and simulation
functions.

# Example

The `joineRmeta` package contains several simulated datasets. One of
these is `simdat` which contains data from 5 studies, each containing
500 simulated individuals. This data first must be transformed into a
joint data object

``` r
library(joineRmeta)
data(simdat)
jointdat<-tojointdata(longitudinal = simdat$longitudinal, 
  survival = simdat$survival, id = "id", longoutcome = "Y", 
  timevarying = c("time","ltime"), survtime = "survtime", cens = "cens",
  time = "time")
jointdat$baseline$study <- as.factor(jointdat$baseline$study)
jointdat$baseline$treat <- as.factor(jointdat$baseline$treat)
```

Using the data held in the `jointdata` object `jointdat` we can then fit
a one stage multi-study joint model.

- `data`: the data object we created above
- `long.formula`: the linear mixed effects model formula for the
  longitudinal sub-model
- `long.rand.ind`: vector of character strings denoting variables to
  assign individual level random effects to
- `long.rand.stud`: vector of character strings denoting variables to
  assign study level random effects to
- `sharingstrct`: the association structure linking the sub-models -
  currently this must be set to `"randprop"`.
- `surv.formula`: the survival formula the survival sub-model.
- `study.name`: character string denoting the variable name of the study
  membership variable.
- `strat`: logical value denoting whether or not the baseline hazard is
  to be stratified by study.

``` r
fit <- jointmeta1(data = jointdat, 
                  long.formula = Y ~ 1 + time + treat, 
                  long.rand.ind = c("int", "time"), 
                  long.rand.stud = "treat",
                  sharingstrct = "randprop",
                  surv.formula = Surv(survtime, cens) ~ treat, 
                  study.name = "study", 
                  strat = T)
```

``` r
summary(fit)
```

Full details on the data and the functions are provided in the help
documentation and package vignette. The purpose of this code is to
simply illustrate the model fitting.

# Funding

This work was supported by the Health eResearch Centre (HeRC) funded by
the Medical Research Council Grant MR/K006665/1.

# Using the latest developmental version

To install the latest **developmental version**, you will need the R
package `devtools` and to run the following code

``` r
library('devtools')
install_github('mesudell/joineRmeta', build_vignettes = FALSE)
```

# References

1.  Cox DR. Regression models and life-tables. *J R Stat Soc Ser B Stat
    Methodol.* 1972; **34(2)**: 187-220.

2.  Gould LA, Boye M, Bois F, et al. Joint modeling of survival and
    longitudinal non-survival data: current methods and issues. Report
    of the DIA Bayesian joint modeling working group. *Stat Med.* 2015;
    **34(14)**: 2181-2195.

3.  Henderson R, Diggle PJ, Dobson A. Joint modelling of longitudinal
    measurements and event time data. *Biostatistics.* 2000; **1(4)**:
    465-480.

4.  Laird NM, Ware JH. Random-effects models for longitudinal data.
    *Biometrics.* 1982; **38(4)**: 963-974.

5.  Wulfsohn MS, Tsiatis AA. A joint model for survival and longitudinal
    data measured with error. *Biometrics.* 1997; **53(1)**: 330-339.
