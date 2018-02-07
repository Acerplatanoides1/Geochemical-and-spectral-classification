# Geochemical and spectral lake core classification 
Developing a method for rapid lake core strata geochemical classification using VIS-NIR Spectroscopy reflectance and geochemical data.

1. Data Preprocessing 

First steps include quantile transformation for Gochemical lake sediment elements as their distribution is non-normal (skewed to left). Raw spectroscopic data is transformed with Savitzky-Golay 1st-derivative although it points out the local maximum peak of the spectrum which is most informative part (represents the elemental overtones- see https://en.wikipedia.org/wiki/Overtone_band). Unfortunetelly transformed spectra have to be additionally aligned in respect to the columns introduced by SG1d filter.




