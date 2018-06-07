Geochemical and spectral reflectance clustering matching for classification of lake core in-depth facies

Introduction

Classification of lake cores and their relation to a changing environment is a matter of Limnological research all over the world. In-depth cores are usually split on a separate layers and analyzed in laboratory to achieve the chemical composition data. Chemical data prepared through standard laboratory measurements eg. spectrometry, wet chemistry and is very time, cost and resource consuming. In other hand spectroscopic measurements using devices such as ASD Inc. VIS-NIR Spectrometer cost less efforts to gather necessary data that could explain the chemical composition of given material (like lake core sediments in this case). Combining both chemical and spectral information about the lake core sediments could be valuable method for a rapid core strata classification not only for Limnological purposes eg. Geology, Paleogeography. Code presented in this repository is an exemplification of classification of chemical and spectral data from in-depth lake core sediments using dimensionallity reduction method Principal Component Analysis (PCA), Gaussian Mixture Modelling and V- Measure model validation method.   

1. Data transformation and spectroscopic data error assessment

The first step to perform the modelling is to transform the data. Data transformation improves the modelling results and is a standard preprocessing step and beforehand PCA chemical data were transformed with quantile normalization and spectral data with Savitzky-Golay first derivative filter (SG1D). Unfortunetelly the SG1D method implemented in package 'Prospectr' introduces additional columns with 0 values to dataset according to derivative used. To surpass that function 'align_lines' were written. After the transformation mean error of spectroscopic data was derived for different window sizes of SG1D transformation and visualized on a graph. Eventually, the best result were taken to procced with models creation - in our case SG1D with window size 7 and 9.  

2. PCA analysis







