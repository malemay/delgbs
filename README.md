# delgbs

## What is delgbs?

*delgbs* is an `R` package which provides tools for detecting copy number variation from genotyping-by-sequencing (GBS) data. *delgbs* bases its CNV calls on the number of reads per sample in discrete bins (e.g. 1-kb bins) located along a reference genome. Given a set of samples, *delgbs* tallies the number of reads per sample and per bin, and then compares the number of reads found in a sample to the population mean by computing a log2 ratio of the relative number of reads in a given bin for a given sample. The log2 ratio profiles thus generated are then segmented using a segmentation algorithm implemented by the package [`robseg`](https://github.com/guillemr/robust-fpop) and segments are categorized as homozygous deletions, hemizygous deletions or duplications according to their average log2 ratio. High-level plotting functions are also provided by the package as a convenience for users.

## Installation

*delgbs* can be installed directly in `R` by calling `devtools::install_github("malemay/delgbs")`. This will directly fetch the package from GitHub and install it on your computer. This requires  [`devtools`](https://cran.r-project.org/web/packages/devtools/index.html) to be installed on your computer. This package is available from CRAN through the usual `ìnstall.packages()` interface.

*delgbs* requires a few [Bioconductor](https://www.bioconductor.org/) packages to be installed on your computer. You can install them by running the following commands in `R`:

```r
{
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite(c("Rsamtools", "S4Vectors", "GenomicRanges", "IRanges", "BiocGenerics"))
}
```

*delgbs* also requires the package *robseg* which can be installed from GitHub by running the command `devtools::install_github("guillemr/robust-fpop")` in `R`.

## Typical usage
The documentation for *delgbs* has yet to be completed. A vignette describing typical usage of the package and example datasets will be uploaded shortly.

## Known issues
None at the moment.

## Citation

If you use this software, please cite:

Lemay, MA., Torkamaneh, D., Rigaill, G. et al. Screening populations for copy number variation using genotyping-by-sequencing: a proof of concept using soybean fast neutron mutants. *BMC Genomics* 20, 634 (2019). [doi:10.1186/s12864-019-5998-1](https://doi.org/10.1186/s12864-019-5998-1)

## Notes
This software is provided without any guarantee. 

Issues, bug reports and questions can be shared on the GitHub page of the project or addressed to the package maintainer (see the package DESCRIPTION for contact information).

