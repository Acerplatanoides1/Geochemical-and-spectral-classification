Geochemical and spectral reflectance clustering matching for classification of lake core in-depth facies


Introduction

Classification of lake cores and their relation to a changing environment is a matter of Limnological research all over the world. In-depth cores are usually split on a separate layers and analyzed in laboratory to achieve the chemical composition data. Chemical data prepared through standard laboratory measurements eg. spectrometry, wet chemistry and is very time, cost and resource consuming. In other hand spectroscopic measurements using devices such as ASD Inc. VIS-NIR Spectrometer cost less efforts to gather necessary data that could explain the chemical composition of given material (like lake core sediments in this case). Combining both chemical and spectral information about the lake core sediments could be valuable method for a rapid core strata classification not only for Limnological purposes eg. Geology, Paleogeography. Code presented in this repository is an exemplification of classification of chemical and spectral data from in-depth lake core sediments using dimensionallity reduction method Principal Component Analysis (PCA), Gaussian Mixture Modelling and V- Measure model validation method.   

1. Data transformation and spectroscopic data error assessment

The first step to perform the modelling is to transform the data. Data transformation improves the modelling results and is a standard preprocessing step and beforehand PCA chemical data were transformed with quantile normalization and spectral data with Savitzky-Golay first derivative filter (SG1D). Unfortunetelly the SG1D method implemented in package 'Prospectr' introduces additional columns with 0 values to dataset according to derivative used. To surpass that function 'align_lines' were written. After the transformation mean error of spectroscopic data was derived for different window sizes of SG1D transformation and visualized on a graph. Eventually, the best result were taken to procced with models creation - in our case SG1D with window size 7 and 9.  

2. PCA analysis

The purpose of PCA analysis is to capture the most of variation of multidimensional dataset and thus reduce the dimensionallity of a given dataset. Chemical data, in our case, contains up to 12 variables (Fe, Mg, K, Na, Ca, Al, Mn, TIC, TOC, SiO2, clay and silt content) and spectral data consist of 2151 variables that are a sequence of wavelengths with it's reflectance values. Usually PCA analysis leads to encapsulation of most important information in max 3 or 4 dimensions, as it is in case of our study dataset. As we possess the matrix with reduced dimensionality, we can procced to Gaussian Mixture Modelling which is a model-based method for clustering of data. 

3. Choosing the distance metrics for clustering

Before our dataset can be modelled we have to choose the distance metrics to be evaluated by a clustering algorithm. With not getting into much of details, we have to transform our data to a distance matrix, which contains distance measures from every observation in given dataset according to choosen method of distance measurement (eg. Euclidean, Manhattan). It is a usual practice to perform calculations for different distance measures and then compare the results, so for our example dataset we choose Manahattan and Mahalanobis distance. As the calculation is done matrices are ready for modelling.

4. Gaussian Mixture Modelling









