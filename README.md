# siteOptiSample 

`siteOptiSample` is an R package developed for the [GEO-TREES](https://geo-trees.org/) initiative to design optimized plot sampling strategies for forest biomass reference sites.

Accurate biomass estimation requires carefully planned ground-truth data collection. Random sampling often misses key structural variations within a forest. `siteOptiSample` generates sampling designs that capture forest heterogeneity within a site. 

`siteOptiSample` can make use of both raster datasets (e.g., airborne LiDAR metrics) or tabular data (e.g., existing tree inventories) to select representative plot locations using a variety of algorithms.

Key use cases include:

* Identifying optimal locations for new tree-inventory plots.
* Selecting representative existing plots for terrestrial LiDAR scanning.

See the [vignette](doc/siteOptiSample.html) for examples of how to use the package.

## Installation

You can install the development version of `siteOptiSample` from GitHub with:

```r
# install.packages("remotes")
remotes::install_github("GEO-TREES/siteOptiSample")
```

## Dependencies

This package relies heavily on the modern R spatial ecosystem. Ensure you have the following installed:

* `terra` - for spatial raster data
* `sf` - for spatial vector data
* `ggplot2` - for visualisation

## License

This project is licensed under the MIT License - see the LICENSE.md file for details.
