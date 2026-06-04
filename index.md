# survFMM

The {survFMM} package fits adaptations of the finite mixture model for
time-to-event endpoints.

Finite mixture models were developed to model continuous or discrete
values arising from separate latent distributions. We have implemented
an adaptation of this modeling approach for time-to-event endpoints,
which require additional consideration due to observations that are
censored and therefore do not have an observed event time.

**IPCW-FMM:** One adaptation is a finite mixture of continuous outcomes,
fit among individuals with an event and correcting for censoring bias
via supplied inverse probability of censoring (IPCW) weights. This
approach and its properties are described in .

**AFT-FMM** The second adaptation is a mixture of accelerated failure
time models (AFT-FMM), in which censored observations are directly
incorporated into the subgroup-specific outcome models. This approach
has more desirable properties when censoring rates are higher, given
that the censoring process does not have to be explicitly modeled and
censored observations are directly incorporated into the log-likelihood.
A manuscript detailing this approach is forthcoming.

## Installation

Install {survFMM} from [GitHub](https://github.com/jalavery/survFMM)
with:

``` r

remotes::install_github("jalavery/survFMM")
```

## References

Jessica A Lavery, Yuan Chen, Katherine S Panageas, Yuanjia Wang,
Unveiling non–small cell lung cancer treatment effect heterogeneity: a
comparative analysis of statistical methods, JNCI: Journal of the
National Cancer Institute, Volume 117, Issue 10, October 2025, Pages
2062–2072, <https://doi.org/10.1093/jnci/djaf176>

Geoffrey J. McLachlan, Sharon X. Lee, Suren I. Rathnayake. 2019. Finite
Mixture Models. Annual Review Statistics and Its Application. 6:355-378.
<https://doi.org/10.1146/annurev-statistics-031017-100325>
