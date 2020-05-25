# dewlap: Comparing dewlap coloration across islands in Caribbean lizards

This R packages accompanies our study on dewlap coloration across habitats and islands in the Caribbean lizard species *Anolis sagrei*. We provide here our data, scripts, results and manuscript. 

## Replication

The analyes of this study were conducted in R, using the RStudio IDE. You can reproduce the study in two ways.

### Option 1

1. Make sure you have the devtools package installed in R, and use `devtools::install_github("rscherrer/dewlap")`. This will install the package of this study and all of the dependencies we used to perform our analyses.
2. Go through the scripts in `analyses` to replicate the study.

### Option 2

1. Download this repository, either the compressed version from the main page or through `git clone`
2. Make sure all the R packages mentioned in the "Imports" and "Suggests" field of the "DESCRIPTION" file are installed on your machine. If not, install them. Packages in "Imports" are used in the functions of this study's package. Packages in "Suggests" are not used in the functions themselves but in the analytical scripts.*
3. Open the project "dewlap.Rproj" in RStudio
4. Under the "Build" tab, click "Install and Restart". This will install the `dewlap` package.
5. You can now go in order through the scripts in `analyses` and replicate the study.

* Note: packages appearing in "Remotes" are not available from CRAN but from GitHub. To install them, make sure you have devtools installed, and use `devtools::install_github("account/package-name")` to install them (e.g. `devtools/install_github("rscherrer/nmgc")`).

## Content

This repository has the default structure of an R package created with RStudio, where `R/` and `man/` contain the package's functions and documentation files. In `data/` you will find our raw data files. In `maps/` we stored the maps we display in the manuscript and scripts used to create them. In `analyses/` are stored all the scripts and outputs (i.e. tables and figures) of our study. The `ms/` folder contains the files related to the manuscript. 

## Manuscript

The manuscript was written in LaTeX using Overleaf and TeXstudio. Note that the TeX files in `ms/` import tables and figures from the `analyses/` folder and is thus not self-contained. This allows the manuscript's tables and figures to automatically update if they are changed from R, for example. The rendered manuscript is available as `main.pdf`. 

### Reviewing (for co-authors)

In order to make use of word processors' reviewing system we also rendered the manuscript as `.docx` and `.odt` documents from `tex` files using [pandoc](https://pandoc.org). So, as a co-author, please use these for reviews, comments and suggestions. The changes will then be incorporated into the original `tex` files. Because file conversion from `tex` to word processor documents is not perfect (in terms of figure resolution and table layout), please also refer to the `pdf` document for full-resolution figures and correct table layout while reviewing text the `docx`/`odt` document. Also please focus on text review only when editing the word processor documents.

