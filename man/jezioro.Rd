\name{jezioro}
\alias{jezioro}
\alias{jezioro}
\docType{package}
\title{
jezioro
}
\description{
The \env{jezioro} package is a collection of functions and datasets intended for use in paleolimnological analyses by members of the Paleoecological Environmental Assessment and Research Laboratory (PEARL).
}

\details{
\tabular{ll}{
Package: \tab jezioro\cr
Type: \tab Package\cr
Version: \tab 0.20\cr
Date: \tab 2020-09-04\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab no\cr
}

\env{jezioro} contains the following functions

\itemize{
         \item \env{cladCount} determines the minimum number of individuals (and their relative abundances) represented in raw cladoceran subfossil count data
         \item \env{dipteranVWHO} harmonize appropriately formatted chironomid count data with the calibration set described in Quinlan and Smol (2001, 2010), construct the described VWHO model, infer VWHO from the chironomid sample data using the model, and perform an analog matching analysis between the sample data and the calibration set 
         \item \env{interpDates} interpolates dates for undated intervals of a dated sediment core
	 \item \env{quickGAM} fit a GAM and perform the related analyses described in Simpson (2018)
         \item \env{vrsChla} infers sediment chlorophyll a concentrations from measurements of visible reflectance
         \item \env{binfordRho} calculates percent water and dry sediment density for individual freeze-dried intervals from a sediment core (values are interpolated for intervals that were not freeze-dried)
         \item \env{binfordActivity} determines corrected 210Pb, 214Bi and 137Cs activities from the appropriate Region of Interest (ROI) values (i.e. output from the gamma counter) and generates a plot of the activities vs. sediment depth
         \item \env{binfordDates} uses the output of \env{binfordRho} and \env{binfordActivity} to calculate 210Pb dates for a sediment core using the constant rate of supply (CRS) model
         \item \env{wiltseBrien} performs a Brien's test on a correlation matrix to test whether changes in the underlying variables are coherent
         \item \env{wiltseDecompose} uses the Brien's test iteratively to identify homogenous and coherent subsets within a correlation matrix
         }

NOTE:\cr
\env{binfordRho}, \env{binfordActivity} and \env{binfordDates} are intended to be run in sequence, as the output from \env{binfordRho} and \env{binfordActivity} provide the input for \env{binfordDates}:
         
\env{jezioro} contains the following data sets:  

\itemize{
         \item \env{chirInfo} a data frame summarizing ecological informatoin for individual chironomid taxa
n example input file for \env{cladCount} (rows are subfossils, columns are samples)
         \item \env{cladCountInput} an example input file for \env{cladCount} (rows are subfossils, columns are samples)
         \item \env{dipteranVWHOInput} an example input file for \env{dipteranVWHO}
         \item \env{vrsChlaInput} an example input file for \env{vrsChla} with sediment intervals as rows and spectrophotometer measurements for wavelengths 650-700 nm, measured every 2 nm
         \item \env{binfordRhoInput} an example input file for \env{binfordRho} formatted as if sediment intervals  were subsampled prior to freeze-drying (i.e. into scintillation vials)
         \item \env{binfordActivityInput} an appropriately formatted input file for \env{binfordActivity}
         \item \env{wiltseInput} an example input file for \env{wiltseBrien} and \env{wiltseDecompose}  (rows are sediment intervals, columns are PCA axis-1 scores for the dominant diatom species in 8 lakes)
         \item \env{tpHall1996} the diatom-TP calibration set described in Hall and Smol (1996)
         \item \env{vwhoQuinlan2010} the dipteran-VWHO calibration set described in Quinlan and Smol (2001, 2010)    
         }
}

\author{
Adam Jeziorski, Joshua Thienpont, Andrew Labaj, Brendan Wiltse

Maintainer: Adam Jeziorski <ajezzz@gmail.com>
}
\references{
Hall RI, Smol JP (1996) Paleolimnological assessment of long-term water-quality changes in south-central Ontario lakes affected by cottage development and acidification. Canadian Journal of Fisheries and Aquatic Sciences 53: 1-17

Simpson GL (2018) Modelling Palaeoecological Time Series Using Generalised Additive Models. Frontiers in Ecology and Evolution DOI:10.3389/fevo.2018.00149

Quinlan R, Smol JP (2001) Chironomid-based inference models for estimating end-of-summer hypolimnetic oxygen from south-central Ontario shield lakes. Freshwater Biology 46: 1529-1551

Quinlan R, Smol JP (2010) Use of Chaoborus subfossil mandibles in models for inferring past hypolimnetic oxygen. Journal of Paleolimnology 44: 43-50.
}

\keyword{ package }

\examples{
#load example cladoceran count data and the calculate number of individuals and relative 
#abundance of each taxon present
data(cladCountInput)
cladCount(cladCountInput)

#load example vrsChla measurements and infer chlorophyll a values
data(vrsChlaInput)
vrsChla(vrsChlaInput)

#load binfordRhoInput data and calculate the percent water of the sediment samples freeze
#-dried in "Whirlpak" bags, while defining custom values for bag weight and core diameter
data(binfordRhoInput)
binfordRho(binfordRhoInput, bagwt=3.75, subsample=FALSE, coreD=3.81)
}
