# Geochemical and spectral reflectance clustering matching for classification of lake core in-depth facies


## Introduction

Classification of lake cores and their relation to a changing environment is a matter of Limnological research all over the world. In-depth cores are usually split on a separate layers and analyzed in laboratory to achieve the chemical composition data. Chemical data prepared through standard laboratory measurements eg. spectrometry, wet chemistry and is very time, cost and resource consuming. In other hand spectroscopic measurements using devices such as ASD Inc. VIS-NIR Spectrometer cost less efforts to gather necessary data that could explain the chemical composition of given material (like lake core sediments in this case). Combining both chemical and spectral information about the lake core sediments could be valuable method for a rapid core strata classification not only for Limnological purposes eg. Geology, Paleogeography. Code presented in this repository is an exemplification of classification of chemical and spectral data from in-depth lake core sediments using dimensionallity reduction method Principal Component Analysis (PCA), Gaussian Mixture Modelling and V- Measure model validation method.   

## 1. Data transformation and spectroscopic data error assessment

The first step to perform the modelling is to transform the data. Data transformation improves the modelling results and is a standard preprocessing step and beforehand PCA chemical data were transformed with quantile normalization and spectral data with Savitzky-Golay first derivative filter (SG1D). Unfortunetelly the SG1D method implemented in package 'Prospectr' introduces additional columns with 0 values to dataset according to derivative used. To surpass that function 'align_lines' were written. After the transformation mean error of spectroscopic data was derived for different window sizes of SG1D transformation and visualized on a graph. Eventually, the best result were taken to procced with models creation - in our case SG1D with window size 7 and 9.  

## 2. PCA analysis

The purpose of PCA analysis is to capture the most of variation of multidimensional dataset and thus reduce the dimensionallity of a given dataset. Chemical data, in our case, contains up to 12 variables (Fe, Mg, K, Na, Ca, Al, Mn, TIC, TOC, SiO2, clay and silt content- all in mg/kg) and spectral data consist of 2151 variables that are a sequence of wavelengths with it's reflectance values. Usually PCA analysis leads to encapsulation of most important information in max 3 or 4 dimensions. In case of our study we will compute the PCA for 70,80 and 90 % of a variation (7,8,9 components) of both chemical and spectral dataset. As we possess the matrix with reduced dimensionality, we can procced to Gaussian Mixture Modelling which is a model-based method for clustering of data. 

## 3. Choosing the distance metrics for clustering

Before our dataset can be modelled we have to choose the distance metrics that will be evaluated by a clustering algorithm. With not getting into much of details, we have to transform our data to a distance matrix, which contains distance measures from every observation in given dataset according to choosen method of distance measurement (eg. Euclidean, Manhattan). It is a usual practice to perform calculations for different distance measures and then compare the results, so for our example dataset we choose Manahattan and Mahalanobis distance. As the calculation is done matrices are ready for modelling.

## 4. Gaussian Mixture Modelling

For chemical and spectral dataset we perform the model-based clustering called Gaussian Mixture Modelling. It is a probabilistic method which bases on Bayes Theorem. Every observation of distance matrix is evaluated prior to assess their relationship to each other (in a distance value). Variables are assumed to have a Gaussian distribution and the algorithm estimates the probability of a given observation that it lies in the certain group (depending on how much groups is choosen). The algorithm can compute the up to 100 different clustering results, but 10 is a sufficient for our goal because after 50 solutions algorithm shows no significant difference in performance. Once we achieve the clustering resultswe can procced to validation.

## 5. Validation of clustering

Validation of clustering is not an easy task and plenty of different methods were developed for this purpose. For our case study  V- Measure have been chosen. The source code was forked from Jesse Mu repository https://gist.github.com/jayelm/9a024359b15516a3d530023ad9ab413d. The chemical and spectral datasets both with Manhattan and Mahalanobis distance and SG1D window size 7 and 9 were compared using this method. The graph in our study case shows that Mahalanobis distance with window size 7 and every amount of PCA components (70,80 90 % of variation of dataset) have the best results 1 that equals the best fit between clustering solutions (1:1).  




