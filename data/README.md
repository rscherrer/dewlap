Most of the data files here are automatically generated from the raw data. Do not edit by hand.

## Content

* `raw`: the folder containing all the raw reflectance spectra
* `counts.csv`: the sample sizes in each habitat on each island
* `island_coordinates.csv`: used to generate the maps
* `reflectance.csv`: the processed and curated reflectance profiles used in the analyses
* `sites.csv`: summarized phenotype (in PC space) for each site on each island
* `specs.rds`: raw reflectance spectra already loaded into R (can be used to save time when replicating the study)

## Dryad

These data were archived on [Dryad](https://datadryad.org). We also uploaded to Dryad the scripts used to process the data from the raw reflectance files and also the scripts to analyze the data and produce the results. All the scripts and data were written within a broader repository, which can be found at [https://github.com/rscherrer/dewlap](https://github.com/rscherrer/dewlap). Hence, if those scripts are to be reused outside of this repository, a prerequisite is that the code is run from the parent directory of the `data` folder. In addition, the scripts automatically save all the resulting figures and tables in various sub-folders of a common `results` folder, which should be on the same level as `data`. Simply running the scripts outside of their original repository, that is, without the right paths and the right folder structure in place, may not work.
