# spida

'''Superseded by [spida2](http://github.com/gmonette/spida2)'''

R package collection of tools for hierarchical and longitudinal data analysis. The package was originally developed to support the Summer Programme in Data Analysis offered at York University from 2000 to 2012. 

The package is being revised and changes are available in the 'spidanew' package. Until 'spidanew' is merged into 'spida' it is
recommended that you install both and load 'spidanew' after 'spida'.

## Installation

To install this package, use the following commands in R:

    if (!require(devtools)) install.packages("devtools")
    library(devtools)
    install_github("gmonette/spida")
    install_github("gmonette/spidanew")

This installs the package from the source, so you will need to have 
R Tools installed on your system.  [R Tools for Windows](https://cran.r-project.org/bin/windows/Rtools/)
takes you to the download page for Windows.  [R Tools for Mac OS X](https://cran.r-project.org/bin/macosx/tools/)
has the required programs for Mac OS X.

To use the package, load it before loading 'spidanew' so the latter will have precedence over 'spida'.

    library(spida)
    library(spidanew)
