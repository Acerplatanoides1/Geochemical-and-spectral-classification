# Geochemical and spectral lake core classification 
Developing a method for rapid lake core strata geochemical classification using VIS-NIR Spectroscopy reflectance and geochemical data.

---
title: "Mixture Spectral Modelling"
output:
  html_notebook: 
    theme: spacelab
  pdf_document: default
---

###########################################Data Transformation#######################################
1. Data setup

```{r setup}
# Zmiana ścieżki do danych
knitr::opts_knit$set(root.dir = normalizePath("C:/Users/Mateusz Obst/Desktop/Hyperspectral RS/R do pracy/Chemia Modelowanie/Modelowanie geochemiczne i spektr/Data/Input/Do analizy/")) 
```

```{r }
# Packages
library(prospectr)
library(ggplot2)
library(mixOmics)
# library(dplyr)
library(mclust)
library(proxy)
library(VGAM)
# library(caret)

# Spectral and Geochemical Data 
g = 'chemia_leb.txt'
s = 'raw_leb.txt'
# load chemical dataset
chemia = as.data.frame(read.table(g, header = TRUE, stringsAsFactors = FALSE))
chemia = chemia[,-1]
#yeo johnson  normalizaion
# a=yeo.johnson(chemia,lambda = 1)
# b=preProcess(chemia,method = "YeoJohnson")


# Quantile normalization for Geochemical data
# source('http://bioconductor.org/biocLite.R')
# biocLite('preprocessCore')
library(preprocessCore)
w = as.matrix(chemia)
d = normalize.quantiles(w) ### works the best
colnames(d) = colnames(chemia)
# load spectra
spektra = as.data.frame(read.table("raw_sprawdzone.txt", header = TRUE, stringsAsFactors = FALSE))
spektra = spektra[,-1]

data = list(spektra = as.matrix(spektra), chemia = as.matrix(d))

```

Aligning function for transformed spectra

```{r}
align_lines <- function(matrix,length,fill=NA) {
  # by rows
  range <- (length-ncol(lines))/2
  m <- matrix(0,nrow(matrix),length)
  m[,range+1:ncol(matrix)+range] <- matrix
  m
}

# szukanie ekstremów spektr
search_for_extrem <- function(x,range) {
  n <- length(x)
  output <- rep(0,n)
  for(i in 1:n) {
    if(is.na(x[i])) next
    left_bound <- ifelse(i-range<1,1,i-range)
    right_bound <- ifelse(i+range>n,n,i+range)
    slice <- x[left_bound:right_bound]
    min_index <- which.min(slice)
    max_index <- which.max(slice)
    if(min_index==range+1 || max_index==range+1) {
      output[i] <- max(slice,na.rm=TRUE)-min(slice,na.rm=TRUE)
      output[i] <- output[i]*ifelse(min_index==range+1,-1,1)
    }
  }
  output
}
```

Spectral SG1D transformation and rows aligning

```{r}

# Transformacja danych spektralnych sg1d
spectra_transform <-function(spectra){
  m = 1
  p = 2
  spectra_sg1d = list(sg1d_w3=savitzkyGolay(spectra,m,p,3),
                      sg1d_w5=savitzkyGolay(spectra,m,p,5),
                      sg1d_w7=savitzkyGolay(spectra,m,p,7),
                      sg1d_w9=savitzkyGolay(spectra,m,p,9),
                      sg1d_w11=savitzkyGolay(spectra,m,p,11),
                      sg1d_w13=savitzkyGolay(spectra,m,p,13),
                      sg1d_w15=savitzkyGolay(spectra,m,p,15),
                      sg1d_w17=savitzkyGolay(spectra,m,p,17),
                      sg1d_w19=savitzkyGolay(spectra,m,p,19),
                      sg1d_w21=savitzkyGolay(spectra,m,p,21),
                      sg1d_w23=savitzkyGolay(spectra,m,p,23),
                      sg1d_w25=savitzkyGolay(spectra,m,p,25),
                      sg1d_w27=savitzkyGolay(spectra,m,p,27),
                      sg1d_w29=savitzkyGolay(spectra,m,p,29),
                      sg1d_w31=savitzkyGolay(spectra,m,p,31))
  spec_ali = list(sg1d_w3 = align_lines(spectra_sg1d[[1]], 2151, fill = NA),
                  sg1d_w5 =align_lines(spectra_sg1d[[2]], 2151, fill = NA),
                  sg1d_w7=align_lines(spectra_sg1d[[3]], 2151, fill = NA),
                  sg1d_w9=align_lines(spectra_sg1d[[4]], 2151, fill = NA),
                  sg1d_w11=align_lines(spectra_sg1d[[5]], 2151, fill = NA),
                  sg1d_w13=align_lines(spectra_sg1d[[6]], 2151, fill = NA),
                  sg1d_w15=align_lines(spectra_sg1d[[7]], 2151, fill = NA),
                  sg1d_w17=align_lines(spectra_sg1d[[8]], 2151, fill = NA),
                  sg1d_w19=align_lines(spectra_sg1d[[9]], 2151, fill = NA),
                  sg1d_w21=align_lines(spectra_sg1d[[10]], 2151, fill = NA),
                  sg1d_w23=align_lines(spectra_sg1d[[11]], 2151, fill = NA),
                  sg1d_w25=align_lines(spectra_sg1d[[12]], 2151, fill = NA),
                  sg1d_w27=align_lines(spectra_sg1d[[13]], 2151, fill = NA),
                  sg1d_w29=align_lines(spectra_sg1d[[14]], 2151, fill = NA),
                  sg1d_w31=align_lines(spectra_sg1d[[15]], 2151, fill = NA))
  return(spec_ali)
}

spec_sg1d = spectra_transform(data[[1]])


```

########################################Calculater matrix error###########################################

Mean error of spectra

```{r}
mean_error_spectra <- function(align_spectra) {
  matrix_error = list(w3_w5 = align_spectra[[1]] + align_spectra[[2]] * 10000,
  w5_w7 = align_spectra[[2]] + align_spectra[[3]] * 10000,
  w7_w9 = align_spectra[[3]] + align_spectra[[4]] * 10000,
  w9_w11 = align_spectra[[4]] + align_spectra[[5]] * 10000,
  w11_w13 = align_spectra[[5]] + align_spectra[[6]] * 10000,
  w13_w15 = align_spectra[[6]] + align_spectra[[7]] * 10000,
  w15_w17 = align_spectra[[7]] + align_spectra[[8]] * 10000,
  w17_w19 = align_spectra[[8]] + align_spectra[[9]] * 10000,
  w19_w21 = align_spectra[[9]] + align_spectra[[10]] * 10000,
  w21_w23 = align_spectra[[10]] + align_spectra[[11]] * 10000,
  w23_w25 = align_spectra[[11]] + align_spectra[[12]] * 10000,
  w25_w27 = align_spectra[[12]] + align_spectra[[13]] * 10000,
  w27_w29 = align_spectra[[13]] + align_spectra[[14]] * 10000,
  w29_w31 = align_spectra[[14]] + align_spectra[[15]] * 10000)
  mean_err = t(data.frame(list(w3=mean(matrix_error[[1]]),
                               w5=mean(matrix_error[[2]]),
                               w7=mean(matrix_error[[3]]),
                               w9=mean(matrix_error[[4]]),
                               w11=mean(matrix_error[[5]]),
                               w13=mean(matrix_error[[6]]),
                               w15=mean(matrix_error[[7]]),
                               w17=mean(matrix_error[[8]]),
                               w19=mean(matrix_error[[9]]),
                               w21=mean(matrix_error[[10]]),
                               w23=mean(matrix_error[[11]]),
                               w25=mean(matrix_error[[12]]),
                               w27=mean(matrix_error[[13]]),
                               w29=mean(matrix_error[[14]]))))
  return(mean_err)
}

mean_error = mean_error_spectra(spec_sg1d)

```

Plot mean error of spectra vs SG1D window size 

```{r}
column = "Średni_błąd"
colnames(mean_error) = column
okna =  data.frame(Okno = c(3, 5, 7, 9, 11, 13,
                            15, 17, 19, 21,23,25,
                            27,29))
wykres1 = cbind(mean_error,okna)

ggplot(wykres1, aes(x = Okno, y = Średni_błąd))+
  geom_point(size = 3, col = "blue", shape = 21)+
  geom_line(col = "blue")+
  theme(plot.title = element_text(hjust = 0.5, size=rel(1.5), face = "bold"))+
  theme(axis.text.x = element_text(face="bold", size=16),
        axis.text.y = element_text(face="bold", size=16))+
  scale_x_continuous(breaks=c(3,5,7,9,11,13,15,17,19,21,23,25,27,29))+
  scale_y_continuous(breaks = c(0.846,0.847,0.848,0.849,0.850,0.851,0.852,0.853))+
  ggtitle("Średni błąd a wielkość okna transformacji spektralnej SG1d")+
  theme_bw()

```

# Continue modelling with window 7 and 9 - best results 
 
#############################################Gaussian Mixture Modelling##################################

PCA Chemistry and Spectra of FULL dataset

1. PCA modelling and extraction of loadings

```{r}
# PCA data sets
chemia = data[[2]]
chemia = as.matrix(chemia[,-1])
w7_w9 = list(spec_sg1d[[3]],spec_sg1d[[4]]) 

# Save PCA results for chemistry and spectra into list

PCA_chemia <- function(chemia) {
  
  PCA_chemia = list(PCA_chem_7= pca(chemia,ncomp = 7,center = F,scale = F),
                    PCA_chem_8= pca(chemia,ncomp = 8,center = F,scale = F),
                    PCA_chem_9= pca(chemia,ncomp = 9,center = F,scale = F),
                    PCA_chem_10= pca(chemia,ncomp = 10,center = F,scale = F))
  
  return(PCA_chemia)
}  

PCA_spectra <- function(spectra) {
  PCA_spectra = list( PCA_w7_7= pca(spectra[[1]],ncomp = 7,center = F,scale = F),
                      PCA_w7_8= pca(spectra[[1]],ncomp = 8,center = F,scale = F),
                      PCA_w7_9= pca(spectra[[1]],ncomp = 9,center = F,scale = F),
                      PCA_w7_10= pca(spectra[[1]],ncomp = 10,center = F,scale = F),
                      PCA_w9_7= pca(spectra[[2]],ncomp = 7,center = F,scale = F),
                      PCA_w9_8= pca(spectra[[2]],ncomp = 8,center = F,scale = F),
                      PCA_w9_9= pca(spectra[[2]],ncomp = 9,center = F,scale = F),
                      PCA_w9_10= pca(spectra[[2]],ncomp = 10,center = F,scale = F))
  return(PCA_spectra)
}

PCA_chem = PCA_chemia(chemia)
PCA_spec = PCA_spectra(w7_w9)

# Save PCA loadings of chemistry and spectra to list 
PCA_load_chem <- function(PCA_chem) {
  PCA_chem_load= list(as.matrix(PCA_chem$PCA_chem_7$x),
                      as.matrix(PCA_chem$PCA_chem_8$x),
                      as.matrix(PCA_chem$PCA_chem_9$x),
                      as.matrix(PCA_chem$PCA_chem_10$x))
  return(PCA_chem_load)
}  

PCA_load_spec <- function(PCA_spec) {
  PCA_spectra_load = list(as.matrix(PCA_spec$PCA_w7_7$x),
                          as.matrix(PCA_spec$PCA_w7_8$x),
                          as.matrix(PCA_spec$PCA_w7_9$x),
                          as.matrix(PCA_spec$PCA_w7_10$x),
                          as.matrix(PCA_spec$PCA_w9_7$x),
                          as.matrix(PCA_spec$PCA_w9_8$x),
                          as.matrix(PCA_spec$PCA_w9_9$x),
                          as.matrix(PCA_spec$PCA_w9_10$x))
  return(PCA_spectra_load)
}

PCA_load_chem = PCA_load_chem(PCA_chem)
PCA_load_spec = PCA_load_spec(PCA_spec)
```

2. Calculation of Manhattan and Mahalanobis Distance 

Chemistry
```{r}
manh_chemia <- function(PCA_loadings) {
  manhat=list(manh_chemia_7= as.matrix(dist(PCA_load_chem[[1]], method = "Manhattan")),
              manh_chemia_8= as.matrix(dist(PCA_load_chem[[2]], method = "Manhattan")),
              manh_chemia_9= as.matrix(dist(PCA_load_chem[[3]], method = "Manhattan")),
              manh_chemia_10= as.matrix(dist(PCA_load_chem[[4]], method = "Manhattan")))
  return(manhat)            
}

mahal_chemia <- function(PCA_loadings) {  
  mahal = list(mahal_chemia_7= as.matrix(mahalanobis(PCA_load_chem[[1]],center = F,cov = cov(PCA_load_chem[[1]]),tol = 1e-25)), 
               mahal_chemia_8= as.matrix(mahalanobis(PCA_load_chem[[2]],center = F,cov = cov(PCA_load_chem[[2]]),tol = 1e-25)), 
               mahal_chemia_9= as.matrix(mahalanobis(PCA_load_chem[[3]],center = F,cov = cov(PCA_load_chem[[3]]),tol = 1e-25)),
               mahal_chemia_10= as.matrix(mahalanobis(PCA_load_chem[[4]],center = F,cov = cov(PCA_load_chem[[4]]),tol = 1e-25)))            
  return(mahal)
}

manh_chem = manh_chemia(PCA_load_chem)
mahal_chem = mahal_chemia(PCA_load_chem)

```

Spectra
```{r}
manh_spectra <- function(PCA_loadings) {
  manhat=list(manh_spec_w7_7= as.matrix(dist(PCA_load_spec[[1]], method = "Manhattan")),
              manh_spec_w7_8= as.matrix(dist(PCA_load_spec[[2]], method = "Manhattan")),
              manh_spec_w7_9= as.matrix(dist(PCA_load_spec[[3]], method = "Manhattan")),
              manh_spec_w7_10= as.matrix(dist(PCA_load_spec[[4]], method = "Manhattan")),
              manh_spec_w9_7= as.matrix(dist(PCA_load_spec[[5]], method = "Manhattan")),
              manh_spec_w9_8= as.matrix(dist(PCA_load_spec[[6]], method = "Manhattan")),
              manh_spec_w9_9= as.matrix(dist(PCA_load_spec[[7]], method = "Manhattan")),
              manh_spec_w9_10= as.matrix(dist(PCA_load_spec[[8]], method = "Manhattan")))
  return(manhat)            
}

mahal_spectra <- function(PCA_loadings) {  
  mahal = list(mahal_spec_w7_7= as.matrix(mahalanobis(PCA_load_spec[[1]],center = F,cov = cov(PCA_load_spec[[1]]),tol = 1e-25)), 
               mahal_spec_w7_8= as.matrix(mahalanobis(PCA_load_spec[[2]],center = F,cov = cov(PCA_load_spec[[2]]),tol = 1e-25)), 
               mahal_spec_w7_9= as.matrix(mahalanobis(PCA_load_spec[[3]],center = F,cov = cov(PCA_load_spec[[3]]),tol = 1e-25)),
               mahal_spec_w7_10= as.matrix(mahalanobis(PCA_load_spec[[4]],center = F,cov = cov(PCA_load_spec[[4]]),tol = 1e-25)),
               mahal_spec_w9_7= as.matrix(mahalanobis(PCA_load_spec[[5]],center = F,cov = cov(PCA_load_spec[[5]]),tol = 1e-25)),
               mahal_spec_w9_8= as.matrix(mahalanobis(PCA_load_spec[[6]],center = F,cov = cov(PCA_load_spec[[6]]),tol = 1e-25)),
               mahal_spec_w9_9= as.matrix(mahalanobis(PCA_load_spec[[7]],center = F,cov = cov(PCA_load_spec[[7]]),tol = 1e-25)),
               mahal_spec_w9_10= as.matrix(mahalanobis(PCA_load_spec[[8]],center = F,cov = cov(PCA_load_spec[[8]]),tol = 1e-25)))
  return(mahal)
}

manh_spect = manh_spectra(PCA_load_chem)
mahal_spect = mahal_spectra(PCA_load_chem)
 
```


#############################Gaussian Mixture Modelling################################################  

1. Mixture model chemia
Manhattan
```{r}
chem_mod_man <- function(PCA_loadings) {
  chem_mod_man = list(chemia_model_man_7= Mclust(manh_chem[[1]],G = 1:10),
                      chemia_model_man_8= Mclust(manh_chem[[2]],G = 1:10),
                      chemia_model_man_9= Mclust(manh_chem[[3]],G = 1:10),
                      chemia_model_man_10= Mclust(manh_chem[[4]],G = 1:10))
  return(chem_mod_man)
}

chem_mod_man = chem_mod_man(manh_chem)

```

Mahalanobis
```{r}
chem_mod_mahal <- function(PCA_loadings) {
  chem_mod_mahal = list(chemia_model_mahal_7= Mclust(mahal_chem[[1]],G = 1:10),
                        chemia_model_mahal_8= Mclust(mahal_chem[[2]],G = 1:10),
                        chemia_model_mahal_9= Mclust(mahal_chem[[3]],G = 1:10),
                        chemia_model_mahal_10= Mclust(mahal_chem[[4]],G = 1:10))
  return(chem_mod_mahal)
}

chem_mod_mahal = chem_mod_mahal(mahal_chem)


```

2. Mixture model spektra
Manhattan
```{r}
spec_mod_manh <- function(PCA_loadings) {
  spec_mod_manh = list(spec_model_man_w7_7= Mclust(manh_spect[[1]],G = 1:10),
                       spec_model_man_w7_8= Mclust(manh_spect[[2]],G = 1:10),
                       spec_model_man_w7_9= Mclust(manh_spect[[3]],G = 1:10),
                       spec_model_man_w7_10= Mclust(manh_spect[[4]],G = 1:10),
                       spec_model_man_w9_7= Mclust(manh_spect[[5]],G = 1:10),
                       spec_model_man_w9_8= Mclust(manh_spect[[6]],G = 1:10),
                       spec_model_man_w9_9= Mclust(manh_spect[[7]],G = 1:10),
                       spec_model_man_w9_10= Mclust(manh_spect[[8]],G = 1:10))
  return(spec_mod_manh)
}

spec_mod_manh = spec_mod_manh(manh_spect)

```

Mahalanobis
```{r}
spec_mod_mahal <- function(PCA_loadings) {
  spec_mod_mahal = list(spec_model_mahal_w7_7= Mclust(mahal_spect[[1]],G = 1:10),
                        spec_model_mahal_w7_8= Mclust(mahal_spect[[2]],G = 1:10),
                        spec_model_mahal_w7_9= Mclust(mahal_spect[[3]],G = 1:10),
                        spec_model_mahal_w7_10= Mclust(mahal_spect[[4]],G = 1:10),
                        spec_model_mahal_w9_7= Mclust(mahal_spect[[5]],G = 1:10),
                        spec_model_mahal_w9_8= Mclust(mahal_spect[[6]],G = 1:10),
                        spec_model_mahal_w9_9= Mclust(mahal_spect[[7]],G = 1:10),
                        spec_model_mahal_w9_10= Mclust(mahal_spect[[8]],G = 1:10))                            
  return(spec_mod_mahal)
}

spec_mod_mahal = spec_mod_mahal(mahal_spect)
spec_mod_mahal
```

Extract model classification

1. Chemistry Manhattan
```{r}
ext_chem_mod_manh <- function(mixture_model) {
  extra = list(chem_clust_manh_7=chem_mod_man[[1]]$classification,
               chem_clust_manh_8=chem_mod_man[[2]]$classification,
               chem_clust_manh_9=chem_mod_man[[3]]$classification,
               chem_clust_manh_10=chem_mod_man[[4]]$classification)
  data = data.frame(extra)
  return(data)
}

ext_chem_mod_manh = ext_chem_mod_manh(chem_mod_man)
ext_chem_mod_manh
```

2. Chemistry Mahalanobis
```{r}
ext_chem_mod_mahal <- function(mixture_model) {
  extra = list(chem_clust_mahal_7=chem_mod_mahal[[1]]$classification,
               chem_clust_mahal_8=chem_mod_mahal[[2]]$classification,
               chem_clust_mahal_9=chem_mod_mahal[[3]]$classification,
               chem_clust_mahal_10=chem_mod_mahal[[4]]$classification)
  data = data.frame(extra)
  return(data)
}

ext_chem_mod_mahal = ext_chem_mod_mahal(chem_mod_mahal)
ext_chem_mod_mahal
```

3. Spectra Manhattan
```{r}
ext_spec_mod_manh <- function(mixture_model) {
  extra = list(spec_clust_manh_w7_7=spec_mod_manh[[1]]$classification,
               spec_clust_manh_w7_8=spec_mod_manh[[2]]$classification,
               spec_clust_manh_w7_9=spec_mod_manh[[3]]$classification,
               spec_clust_manh_w7_10=spec_mod_manh[[4]]$classification,
               spec_clust_manh_w9_7=spec_mod_manh[[5]]$classification,
               spec_clust_manh_w9_8=spec_mod_manh[[6]]$classification,
               spec_clust_manh_w9_9=spec_mod_manh[[7]]$classification,
               spec_clust_manh_w9_10=spec_mod_manh[[8]]$classification)
  data = data.frame(extra)
  return(data)
}

ext_spec_mod_manh = ext_spec_mod_manh(spec_mod_manh)
ext_spec_mod_manh
```

4. Spectra Mahalanobis
```{r}
ext_spec_mod_mahal <- function(mixture_model) {
  extra = list(spec_clust_mahal_w7_7=spec_mod_mahal[[1]]$classification,
               spec_clust_mahal_w7_8=spec_mod_mahal[[2]]$classification,
               spec_clust_mahal_w7_9=spec_mod_mahal[[3]]$classification,
               spec_clust_mahal_w7_10=spec_mod_mahal[[4]]$classification,
               spec_clust_mahal_w9_7=spec_mod_mahal[[5]]$classification,
               spec_clust_mahal_w9_8=spec_mod_mahal[[6]]$classification,
               spec_clust_mahal_w9_9=spec_mod_mahal[[7]]$classification,
               spec_clust_mahal_w9_10=spec_mod_mahal[[8]]$classification)
  data = data.frame(extra)
  return(data)
}

ext_spec_mod_mahal = ext_spec_mod_mahal(spec_mod_mahal)




```


###############################################Clustering Validation##########################

1. V- Measure validation function - comparison of chemical and spectral classifications

```{r}
library(infotheo)

v.measure <- function(a, b) {
  mi <- mutinformation(a, b)
  entropy.a <- entropy(a)
  entropy.b <- entropy(b)
  if (entropy.a == 0.0) {
    homogeneity <- 1.0
  } else {
    homogeneity <- mi / entropy.a
  }
  if (entropy.b == 0.0) {
    completeness <- 1.0
  } else {
    completeness <- mi / entropy.b
  }
  if (homogeneity + completeness == 0.0) {
    v.measure.score <- 0.0
  } else {
    v.measure.score <- (2.0 * homogeneity * completeness
                       / (homogeneity + completeness))
  }
  # Can also return homogeneity and completeness if wanted
  v.measure.score
}
```

2. Functions for V-Measure calculation

```{r}
# V-Measure Manhattan 
v.measure_manh <- function(clustering_chemia,clustering_spectra) {
  v_measure_chem_spec_manh = t(data.frame(
    v.mea_manh7_w7_7 = v.measure(ext_chem_mod_manh$chem_clust_manh_7,ext_spec_mod_manh$spec_clust_manh_w7_7),
                                 v.mea_manh8_w7_8 =v.measure(ext_chem_mod_manh$chem_clust_manh_8,ext_spec_mod_manh$spec_clust_manh_w7_8),
                                 v.mea_manh9_w7_9 =v.measure(ext_chem_mod_manh$chem_clust_manh_9,ext_spec_mod_manh$spec_clust_manh_w7_9),
                                 v.mea_manh10_w7_10 =v.measure(ext_chem_mod_manh$chem_clust_manh_10,ext_spec_mod_manh$spec_clust_manh_w7_10),
                                 v.mea_manh7_w9_7 =v.measure(ext_chem_mod_manh$chem_clust_manh_7,ext_spec_mod_manh$spec_clust_manh_w9_7),
                                 v.mea_manh8_w9_8 =v.measure(ext_chem_mod_manh$chem_clust_manh_8,ext_spec_mod_manh$spec_clust_manh_w9_8),
                                 v.mea_manh9_w9_9 =v.measure(ext_chem_mod_manh$chem_clust_manh_9,ext_spec_mod_manh$spec_clust_manh_w9_9),
                                 v.mea_manh10_w9_10 =v.measure(ext_chem_mod_manh$chem_clust_manh_10,ext_spec_mod_manh$spec_clust_manh_w9_10)))
    colnames(v_measure_chem_spec_manh) = "V-Measure Manh"
    return(v_measure_chem_spec_manh)
}
 
v.measure_manh =v.measure_manh(ext_chem_mod_manh,ext_spec_mod_manh)
v.measure_manh
```

```{r}
# V-Measure Mahalanobis
v.measure_mahal <- function(clustering_chemia,clustering_spectra) {
  v_measure_chem_spec_mahal = t(data.frame(
    v.mea_mahal7_w7_7 = v.measure(ext_chem_mod_mahal$chem_clust_mahal_7,ext_spec_mod_mahal$spec_clust_mahal_w7_7),
                                 v.mea_mahal8_w7_8 =v.measure(ext_chem_mod_mahal$chem_clust_mahal_8,ext_spec_mod_mahal$spec_clust_mahal_w7_8),
                                 v.mea_mahal9_w7_9 =v.measure(ext_chem_mod_mahal$chem_clust_mahal_9,ext_spec_mod_mahal$spec_clust_mahal_w7_9),
                                 v.mea_mahal10_w7_10 =v.measure(ext_chem_mod_mahal$chem_clust_mahal_10,ext_spec_mod_mahal$spec_clust_mahal_w7_10),
                                 v.mea_mahal7_w9_7 =v.measure(ext_chem_mod_mahal$chem_clust_mahal_7,ext_spec_mod_mahal$spec_clust_mahal_w9_7),
                                 v.mea_mahal8_w9_8 =v.measure(ext_chem_mod_mahal$chem_clust_mahal_8,ext_spec_mod_mahal$spec_clust_mahal_w9_8),
                                 v.mea_mahal9_w9_9 =v.measure(ext_chem_mod_mahal$chem_clust_mahal_9,ext_spec_mod_mahal$spec_clust_mahal_w9_9),
                                 v.mea_mahal_10_w9_10 =v.measure(ext_chem_mod_mahal$chem_clust_mahal_10,ext_spec_mod_mahal$spec_clust_mahal_w9_10)))
    colnames(v_measure_chem_spec_mahal) = "V-Measure Mahal"
    return(v_measure_chem_spec_mahal)
}
?mahalanobis()
v.measure_mahal =v.measure_mahal(ext_chem_mod_mahal,ext_spec_mod_mahal)
v.measure_mahal

```

3. V-Measure graph
```{r}
v.measure_manh
v.measure_mahal

v.measure_a = rbind(v.measure_manh,v.measure_mahal)
e = c("Manhattan7_w7","Manhattan8_w7","Manhattan9_w7","Manhattan10_w7","Manhattan7_w9","Manhattan8_w9","Manhattan9_w9","Manhattan10_w9",
      "Mahalanobis7_w7","Mahalanobis8_w7","Mahalanobis9_w7","Mahalanobis10_w7","Mahalanobis7_w9","Mahalanobis8_w9","Mahalanobis9_w9","Mahalanobis10_w9")

v.measure_b = data.frame(e)


v.measure_score = data.frame(v.measure_a, labels = v.measure_b)
colnames(v.measure_score) = c("V.Measure","Dystans_i_wielkosc_okna")
v.measure_score[,"V.Measure"] = round(v.measure_score[,"V.Measure"],digits =  3)

ggplot(v.measure_score, aes(x = V.Measure,y = Dystans_i_wielkosc_okna))+
  geom_point(size = 4, col = "blue", shape = 20)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size=rel(1.5), face = "bold"))+
  theme(axis.text.x = element_text(face="bold", size=11),
        axis.text.y = element_text(face="bold", size=11))+
  scale_x_continuous(breaks=c(0.4,0.5,0.6,0.7,0.8,0.9,1))+
  geom_text(aes(label=V.Measure,hjust = 0,vjust = 1), size=4, colour="red")+
  ggtitle("V-Measure dla dystansu Manhattańskiego i Mahalanobisa")
  
  
  

```

