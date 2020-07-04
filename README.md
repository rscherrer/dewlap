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

The manuscript was written in LaTeX using Overleaf and TeXstudio. The figures and tables for the manuscript are in the folders `figures`, `suppfigures`, `tables` and `supptables`. Run the `fetch_content.sh` script to copy the relevant figures and tables from the `analyses` folder into the `ms` folder. The rendered manuscript is available as `main.pdf`. The script `run_pandoc.sh` converts the TeX files into two word-processor document versions of the manuscript, `scherrer2020.docx` and `scherrer2020.odt`. The ODT version was uploaded on Google Drive and shared with collaborators, and this is where the review of the draft was mostly done.
