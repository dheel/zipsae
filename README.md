
<!-- README.md is generated from README.Rmd. Please edit that file -->

# zipsae

<!-- badges: start -->
<!-- badges: end -->

This function produces empirical best linier unbiased predictions
(EBLUPs) for Zero-Inflated data and its Relative Standard Error. Small
Area Estimation with Zero-Inflated Model (SAE-ZIP) is a model developed
for Zero-Inflated data that can lead us to overdispersion situation. To
handle this kind of situation, this model is created. The model in this
package is based on Small Area Estimation with Zero-Inflated Poisson
model proposed by Dian Christien Arisona
(2018)<https://repository.ipb.ac.id/handle/123456789/92308>. For the
data sample itself, we use combination method between Roberto Benavent
and Domingo Morales (2015)<doi:10.1016/j.csda.2015.07.013> and Sabine
Krieg, Harm Jan Boonstra and Marc Smeets
(2016)<doi:10.1515/jos-2016-0051>.

## Authors

Fadheel Wisnu Utomo, Ika Yuni Wulansari

## Maintainer

Fadheel Wisnu Utomo <221709671@stis.ac.id>

## Installation

You can install the released version of zipsae from
[CRAN](https://CRAN.R-project.org) or find on my github repository
[Github](https://github.com/dheel)

## Example

``` r
##load the dataset in package
library(zipsae)
data(dataSAEZIP)

##Extract the vardir (sampling error)
dataSAEZIP$vardir -> sError

##Compute the data with SAE ZIP model
formula = (y~x1)
zipsae(data = dataSAEZIP, vardir = sError, formula) -> saezip

head(saezip$estimate)
#>           [,1]
#> [1,] 0.2925708
#> [2,] 0.2790501
#> [3,] 0.2772425
#> [4,] 0.2884874
#> [5,] 0.2931530
#> [6,] 0.2970365
## saezip$estimate        #to see the result of Small Area Estimation with Zero-Inflated Model
## saezip$dispersion$rse  #to see the relative standard error from the estimation
## saezip$coefficient$lambda   #to see the estimator which is gained from the non-zero compilation data.
## saezip$coefficient$omega   #to see the estimator which is gained from the complete compilation data.
```

## References

-   Arisona, D.C. (2018). Kajian Pendugaan Area Kecil pada Data
    Overdispersi Menggunakan Regresi Zero-Inflated Poisson. Bogor: Bogor
    Agricultural University.
-   Benavent, Roberto & Morales, Domingo. (2015). Multivariate
    Fay-Herriot models for small area estimation. Computational
    Statistics and Data Analysis 94 2016 372-390. DOI:
    10.1016/j.csda.2015.07.013.
-   Rao, J.N.K & Molina. (2015). Small Area Estimation 2nd Edition. New
    York: John Wiley and Sons, Inc.
-   S. Krieg, H. J. Boonstra, and M. Smeets. Small-area estimation with
    zero-inflated data – a simulation study. J. Off. Stat., vol. 32, no.
    4, pp. 963–986, 2016, doi: 10.1515/JOS-2016-0051
