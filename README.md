# Geochemical and spectral lake core classification 
Developing a method for rapid lake core strata geochemical classification using VIS-NIR Spectroscopy reflectance and elemental data.

1. Data Preprocessing 

First steps include quantile transformation for Gochemical lake sediment data as their distribution is non-normal (skewed to left or right). Raw spectroscopic data is transformed with Savitzky-Golay 1st-derivative (SG1d), although it points out the local maximum peak of the spectrum which is most informative part (represents the elemental overtones- see https://en.wikipedia.org/wiki/Overtone_band). SG1d performance relies on the filter kernel window size parameter, which always is an odd number (3x3, 5x5, 7x7 etc.). Unfortunetelly, transformed spectra have to be additionally aligned in respect to the additional columns introduced by SG1d filter function. To resolve the right window size we can calculate the mean error of different window settings. This allows us to distinguish proper window size to choose for further analysis. In our case window size 7x7 and 9x9 indicate the best results and will be used in next steppes.          

2. Data Analysis

a) 




