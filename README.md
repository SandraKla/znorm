# *znorm*: Box-Cox standardization of positive laboratory measurements

The *znorm* package provides a Box-Cox standardization of positive laboratory measurements using reference intervals, with the corresponding inverse transformation.

## Installation

You can install `znorm` directly from GitHub. The snippet below automatically checks and installs required dependencies:

```r
if ("znorm" %in% rownames(installed.packages())) {
  library(znorm)
} else {
  if ("devtools" %in% rownames(installed.packages())) {
    library(devtools)
  } else {
    install.packages("devtools")
    library(devtools)
  }
  devtools::install_github("SandraKla/znorm")
  library(znorm)
}
```

Once the installation is complete, load the package into your R session using the following command:
```r
library(znorm)
```

The package will then be ready for use in your R environment (e.g. in RStudio). To see the documentation of the package with all its help files, please enter

```r
help(package = znorm)
```

## Usage

Here is a basic example demonstrating how to use *znorm*:

```r
library(znorm)
albumin <- c(42, 34, 38, 43, 50, 42, 27, 31, 24)
scores <- znorm(albumin, limits = c(35, 52), lambda = 0.5)
scores
iznorm(scores, limits = c(35, 52), lambda = 0.5)

x <- c(albumin = 42, bilirubin = 8)
limits <- rbind(c(35, 52), c(2, 21))
znorm(x, limits, lambda = c(1, 0.33))

```

## Publication

*zlog* has been published in:  Hoffmann et al. (2017), [doi:10.1515/labmed-2016-0087](https://doi.org/10.1515/labmed-2016-0087)