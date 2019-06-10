# jezioro

The jezioro R package is a collection of functions and datasets intended for use by members of the Paleoecological Environmental Assessment and Research Laboratory (PEARL).

The package includes the following functions (along with example data to illustrate their use):
* *cladCount*: determine the number of individuals (and their relative abundances) represented in raw cladoceran subfossil data
* *dipteranVWHO*: infer VWHO from dipteran relative abundances using the dipteran-inferred VWHO calibration set described in Quinlan and Smol (2001; 2010)
* *interpDates*: interpolate dates for the undated intervals of a dated sediment core
* *vrsChla*: infer chlorophyll *a* concentration in lake sediment from visible reflectance measurements
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
All dependencies are related to the various examples in the general guide vignette and include: *dplyr*, *ggplot2*, *knitr*, *magrittr*, *mapproj*, *rgdal*, *rioja*, and *rmarkdown*.


## Installation
The most recent version of the *jezioro* package can be downloaded from the "release" tab. Both source code and compiled binaries (for windows + macOS/linux) are available there as assets. These binaries can then be installed as ["a local file"](https://www.rdocumentation.org/packages/utils/versions/3.5.1/topics/install.packages).

Alternately, *jezioro* can be installed directly into a local R envrionment from the github repository using the *devtools* package. However, note that this method will not provide a local copy of vignettes/guides.
```
library("devtools")
install_github("shiggo/jezioro")
```
