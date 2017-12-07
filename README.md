# Local_in_time_statistics

This repository contains the code used for local in time statistics analyses of gene expression data with respect to time to diagnosis as described in [Holden et al. 2017.](http://www.ingentaconnect.com/content/doaj/11799870/2017/00000055/00000001/art00003)


# Repository

This method requires R and some dependencies mentioned in the scripts.

The five scripts in this repository constitute the analyses performed using the breast cancer case-control data set described in the publication.
The scripts are adopted from the original analysis to analyse any prospective case-control data set.


# Usage

This code is adopted to analyse case-control differences in gene expression levels in cancer studies employing a prospective design.
Data on time to diagnosis must be available and the methods allow for comparisons between two strata (e.g. with and without spread).

The main script is analyseCancerData.r which calls the other scripts. 
The analysis assumes a data file containing microarray-based gene expression data in …data/prospectiveCancerData.R.
Also, adjustments of e.g. reading of the data, variable names and file paths must be performed for each data set.
Further details can be found in each script.


# Citation

When used, the method should be cited with reference to Holden et al. 2017.

Also, the code is to be cited using the following doi: doi.org/10.2147/AGG.S130004 when describing the methods or using e.g. the following sentence as an acknowledgement: "The code used for analysis was provided in a repository at https://github.com/theresenoest/Local_in_time_statistics".

Lisence: MIT


# Issues

Please report any issues regarding the contents of the repository at the issues site.
