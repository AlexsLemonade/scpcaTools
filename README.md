
# scpcaTools: Useful tools for analysis of single-cell RNA seq counts data

  <!-- badges: start -->
  [![R-CMD-check](https://github.com/AlexsLemonade/scpcaTools/workflows/R-CMD-check/badge.svg)](https://github.com/AlexsLemonade/scpcaTools/actions)
  <!-- badges: end --> 

The `scpcaTools` package contains a set of tools for working with single-cell and single-nuclei RNA-seq counts data.
Mainly, this package can work with data that has been produced using Alevin, Alevin-Fry, Cellranger, or Kallisto. 
Counts matrices from either of these pre-processing tools can be imported into R using the `import_quant_data` function to return a `SingleCellExperiment` object to be used for downstream analysis. 

Currently, the `import_quant_data` function can support single-nuclei RNA seq counts data that has been aligned to a pre-mRNA index by using the `intron_mode=TRUE` and `which_counts=intron` options. There is also support to import data processed using Alevin-fry's USA mode by using the `usa_mode=TRUE` option. 

## Installation 

The package can be installed from github with:

```r
remotes::install_github("AlexsLemonade/scpcaTools")
```


## Testing

To test the import functions in this package, you must also have installed the [`scpcaData` package](https://github.com/AlexsLemonade/scpcaData).
This can be done by running `remotes::install_github("AlexsLemonade/scpcaData")`
