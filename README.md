# jezioro

The jezioro R package is a collection of functions and datasets intended for use by members of the Paleoecological Environmental Assessment and Research Laboratory (PEARL).

The package includes the following functions (along with example data to illustrate their use):
* *cladCount* - determines the number of individuals (and their relative abundances) represented in raw cladoceran subfossil data
interpDates - interpolates dates for the undated intervals of a dated sediment core
* *vrsChla* - infers lake sediment chlorophyll a concentrations from measurements of visible reflectance
* the 'binford' functions (*binfordRho*, *binfordActivity*, *binfordDates*) collectively calculate 210Pb dates for a sediment core using the methods described in Binford 1990, Schelske 1994, and Appleby 2001.
* the 'wiltse' functions (*wiltseBrien*, *wiltseDecompose*) test for homogeneity and coherence within a dataset. They are slightly modified (so that they no longer write to the global environment) versions of the ‘brien’ and ‘decompose’ functions described in Brendan Wiltse’s 2014 PhD thesis.

The package also contains some commonly used calibration sets to facilitate the construction of transfer functions for application to individual sediment cores:
* *tpHall1996* is the diatom-inferred TP calibration set of 54 south-central Ontario lakes described in Hall and Smol (1996)
* *vwhoQuinlan2010* is the dipteran-inferred VWHO calibration set of 54 south-central Ontario lakes described in Quinlan and Smol (2001; 2010) 