This folder contains the R scripts used in our study.
To replicate our study, the scripts should be run in the order indicated by their numbers.

## Content

* `01_pocess_data.R`: a script to curate, process and assemble the raw spectral data into a dataset ready for analysis
* `02_run_analyses.R`: a script running all the analyses of our study on the processed data and saving various figures and tables as output
* `03_convert_tables.R`: a script to convert the tables saved by the previous script into LaTeX tables, ready to plug into the manuscript
* `04_median_distance.R`: a script to compute the median distance between localities within islands based on latitudes and longitudes

Please refer to the headers and comments within each script for more details about the procedures implemented.
You will also find a detailed description of our processing steps in the README file inside the `../data` folder.

## Disclaimer

These scripts were archived on [Dryad](https://datadryad.org), separately from the data they use (`../data`). 
Hence, if you found this README on Dryad please be aware that the scripts were written as part of a broader repository (available at [https://github.com/rscherrer/dewlap](https://github.com/rscherrer/dewlap)) with a specific folder structure.
For example, the scripts should be run from a working directory that contains the present `data` directory as well as a `results` folder containing specific sub-folders where the scripts are supposed to save tables and figures to.
Simply running the scripts outside of their original repository may not work, or you may have to update the paths in the scripts to make them work. 
