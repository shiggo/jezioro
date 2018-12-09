# jezioro

The jezioro R package is a collection of functions and datasets intended for use by members of the Paleoecological Environmental Assessment and Research Laboratory (PEARL).

The package includes the following functions (along with example data to illustrate their use):
* *cladCount*: determine the number of individuals (and their relative abundances) represented in raw cladoceran subfossil data
* *interpDates*: interpolate dates for the undated intervals of a dated sediment core
* *vrsChla*: infer chlorophyll *a* concentration in lake sediment from visible reflectance measurments
* the 'binford' functions collectively calculate <sup>210</sup>Pb dates for a sediment core using the methods described in Binford 1990, Schelske 1994, and Appleby 2001
  * *binfordRho*: determine water content and dry sediment density for freeze-dried intervals of a sediment core
  * *binfordActivity*: calculate <sup>210</sup>Pb, <sup>137</sup>Cs, and <sup>214</sup>Bi activities corrected for the efficiency of the gamma counter
  * *binfordDates*: use the unsupported fraction of <sup>210</sup>Pb activity in sediment samples to calculate sediment ages via the Constant Rate of Supply (CRS) model 
* the 'wiltse' functions test for homogeneity and coherence within a dataset. They are slightly modified versions of the 'brien' and 'decompose' functions described in Brendan Wiltse’s 2014 PhD thesis.
  * *wiltseBrien*: perform the Brien test (Brien et al. 1984) on variables contained within a data frame
  * *wiltseDecompose*: build upon the Brien test by identifying homogenous and coherent subsets within the correlation matrix

The package also contains some commonly used calibration sets to facilitate the construction of transfer functions for application to individual sediment cores:
* *tpHall1996*: diatom-inferred TP calibration set of 54 south-central Ontario lakes described in Hall and Smol (1996)
* *vwhoQuinlan2010*: dipteran-inferred VWHO calibration set of 54 south-central Ontario lakes described in Quinlan and Smol (2001; 2010)


## Documentation
The package also contains two guides or vignettes:

* The [user guide](https://shiggo.github.io/jezioro/vignettes/jezioroGuide.html) provides detailed information on the functions and data sets within the package as well as examples of their use.

* A [general guide to the use of R at PEARL](https://shiggo.github.io/jezioro/vignettes/RGuide.html) is also present. This guide is an ongoing attempt to provide an overview of the most commonly used functions and packages relevant to the analyses performed at PEARL.


## Dependencies
The package dependencies are due to running the examples in the overview vignetter and include: *dplyr*, *ggplot2*, *knitr*, *magrittr*, *mapproj*, *rgdal*, *rioja*, and *rmarkdown*.


## Building and Installing
The simplest way to obtain a working install of *jezioro* is directly from the github repository using the *devtools* package for the R sofware environment. However, this method will not provide a local copy of vignettes/guides.
```
library("devtools")
install_github("shiggo/jezioro")
```
Alternately, download the source (green button on the github page), using the "Download ZIP" option, extract the archive and run the build function pointing the pkg argument at the extracted folder to produce a zip binary version of the package that can be installed in your local R environment using the "install from a local archive" option.
