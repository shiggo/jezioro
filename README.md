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
The [jezioro User Guide](https://shiggo.github.io/jezioro/vignettes/jezioroGuide.html) contains detailed information on the package contents and their use.

The package also contains [An Introduction to using R at PEARL](https://shiggo.github.io/jezioro/vignettes/RGuide.html) which is a general overview of how to use R to perform the more common analyses performed at PEARL.


## Dependencies
The package dependencies are due to running the examples in the overview vignetter and include: *dplyr*, *ggplot2*, *knitr*, *magrittr*, *mapproj*, *rgdal*, *rioja*, and *rmarkdown*.


## Compiling and Installing
To build *jezioro*, first download and extract the jezioro-master.zip archive, then use the *build* command provided by the R *devtools* package to build a local package archive to install.