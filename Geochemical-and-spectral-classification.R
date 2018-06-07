# Spectral Mixture Modelling with Geochemical Data

#Data Transformation
1. Data setup

# Load package 
library(prospectr)
library(ggplot2)
library(mixOmics)
library(mclust)
library(proxy)
library(VGAM)

#Quantile normalization
# source('http://bioconductor.org/biocLite.R')
# biocLite('preprocessCore')
library(preprocessCore)

# Load data files
g = 'chemia_leb.txt'
s = 'raw_spectra.txt'

# Load chemical dataset 
chemia = data.frame(read.table(g, header = TRUE, stringsAsFactors = FALSE))
id_chemia = chemia[,1]
chemia = chemia[,-1]

# Normalize dataset
chemia = normalize.quantiles(as.matrix(chemia)) 
chemia = round(chemia, digits = 2)

# Load spectral dataset 
spectra = as.matrix(read.table(s, header = TRUE, stringsAsFactors = FALSE, row.names = 1))
spectra = round(spectra, digits = 7)

#Aligning function for transformed spectra

align_lines <- function(lines,length,fill=NA) {
   # by rows
   range <- (length-ncol(lines))/2
   m <- matrix(0,nrow(lines),length)
   m[,range+1:ncol(lines)+range] <- lines
   m
}

#Spectral SG1D transformation and rows aligning

m=1
p=2

aligned = numeric(0)
for ( i in seq(3,31,2)) {
   tmp <- savitzkyGolay(spectra,m,p,i)
   tmpl <- align_lines(tmp, 2151)
   aligned <- append(aligned,list(tmpl))
}

#Mean error of spectra

err = numeric(0)
error = numeric(0)
for (i in 1:14) {
  er = aligned[[i]] + aligned[[i+1]] * 10000
  err = append(err,list(er))
  error[[i]] = mean(err[[i]])
}

#Plot mean error of spectra vs SG1D window size 

column = "mean_error"
error = data.frame(error)
colnames(error) = column
okna =  data.frame(window = c(seq(3,29,2)))
wykres1 = cbind(error,window)

ggplot(plot1, aes(x = window, y = mean_error))+
  geom_point(size = 3, col = "blue", shape = 21)+
  geom_line(col = "blue")+
  scale_x_continuous(breaks=c(seq(3,29,2)))+
  scale_y_continuous(breaks = c(0.846,0.847,0.848,0.849,0.850,0.851,0.852,0.853))+
  ggtitle("Mean error and window size after transformation of SG1d")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size=rel(1.5), face = "bold"))+
  theme(axis.text.x = element_text(face="bold", size=10),
        axis.text.y = element_text(face="bold", size=10))

#Continue modelling with window 7 and 9 - best results. 
 
#Gaussian Mixture Modelling

PCA Chemistry and Spectra 

1. PCA modelling and extraction of loadings

# PCA for chemistry

PCA_chem = numeric(0)
for (i in 7:10) {
  PCA_che = pca(chemia,ncomp = i,center = F,scale = F)
  PCA_chem =append(PCA_chem,PCA_che) 
  }

# Save PCA loadings for chemistry into list

PCA_chem_load = list(PCA_chem[[10]],PCA_chem[[25]],PCA_chem[[40]],PCA_chem[[55]])

# PCA for spectra

PCA_spectra_w7 = numeric(0)
PCA_spectra_w9 = numeric(0)
for (i in 7:10) {
  PCA_spec_w7 = pca(aligned[[3]],ncomp = i,center = F,scale = F) # for windows size 7
  PCA_spec_w9 = pca(aligned[[4]],ncomp = i,center = F,scale = F) # for windows size 9
  PCA_spectra_w7 =append(PCA_spectra_w7,PCA_spec_w7) 
  PCA_spectra_w9 =append(PCA_spectra_w9,PCA_spec_w9) 
  }
# Save PCA loadings for spectra into list

PCA_spect_load = list(PCA_spectra_w7[[10]],PCA_spectra_w7[[25]],PCA_spectra_w7[[40]],PCA_spectra_w7[[55]],
                      PCA_spectra_w9[[10]],PCA_spectra_w9[[25]],PCA_spectra_w9[[40]],PCA_spectra_w9[[55]])


2. Calculation of Manhattan and Mahalanobis Distance with PCA Loadings

#Chemistry

# Manhattan Distance Chemistry

manh_chemistry = numeric(0)
for (i in 1:4) {
  manh_chem= as.matrix(dist(PCA_chem_load[[i]], method = "Manhattan"))
  manh_chemistry = append(manh_chemistry, list(manh_chem))
}

# Mahalanobis Distance Chemistry

mahal_chemistry = numeric(0)
for (i in 1:4) {
  mahal_chem =  as.matrix(mahalanobis(PCA_chem_load[[i]],center = F,cov = cov(PCA_chem_load[[i]]),tol = 1e-25))
  mahal_chemistry = append(mahal_chemistry, list(mahal_chem))
}

#Spectra

# Manhattan Distance Spectra

manh_spect_w7 = numeric(0)
manh_spect_w9 = numeric(0)
for (i in 1:4) {
  manh_spec= as.matrix(dist(PCA_spect_load[[i]], method = "Manhattan"))
  manh_spect_w7 = append(manh_spect_w7, list(manh_spec))
}

for (i in 5:8) {
  manh_spec= as.matrix(dist(PCA_spect_load[[i]], method = "Manhattan"))
  manh_spect_w9 = append(manh_spect_w9, list(manh_spec))
  }

manh_spect = list(manh_spect_w7[[1]],manh_spect_w7[[2]],manh_spect_w7[[3]],manh_spect_w7[[4]],
                  manh_spect_w9[[1]],manh_spect_w9[[2]],manh_spect_w9[[3]],manh_spect_w9[[4]])

# Mahalanobis Distance Spectra
mahal_spect_w7 = numeric(0)
mahal_spect_w9 = numeric(0)
for (i in 1:4) {
  mahal_spec =  as.matrix(mahalanobis(PCA_spect_load[[i]],center = F,cov = cov(PCA_spect_load[[i]]),tol = 1e-25))
  mahal_spect_w7 = append(mahal_spect_w7, list(mahal_spec))
  }

for (i in 5:8) {
  mahal_spec= as.matrix(dist(PCA_spect_load[[i]], method = "Manhattan"))
  mahal_spect_w9 = append(manh_spect_w9, list(manh_spec))
  }

mahal_spect = list(mahal_spect_w7[[1]],mahal_spect_w7[[2]],mahal_spect_w7[[3]],mahal_spect_w7[[4]],
                  mahal_spect_w9[[1]],mahal_spect_w9[[2]],mahal_spect_w9[[3]],mahal_spect_w9[[4]])

Gaussian Mixture Modelling  

1. Mixture model chemistry

#Manhattan

chem_man_mod = numeric(0)
for (i in 1:4) {
  chem_man_mo= Mclust(manh_chemistry[[i]],G = 1:10)
  chem_man_mod = append(chem_man_mod,chem_man_mo)
  }

#Mahalanobis

chem_mahal_mod = numeric(0)
for (i in 1:4) {
  chem_maha_mo= Mclust(mahal_chemistry[[i]],G = 1:10)
  chem_mahal_mod = append(chem_mahal_mod,chem_maha_mo)
  }

# Clustering results
chem_mod_man = list(chem_man_mod[[14]],chem_man_mod[[29]],chem_man_mod[[44]],chem_man_mod[[59]])
                
chem_mod_mahal = list(chem_mahal_mod[[14]],chem_mahal_mod[[29]],chem_mahal_mod[[44]],chem_mahal_mod[[59]])

2. Mixture model spectra

#Manhattan

spect_man_mod_w7 = numeric(0)
for (i in 1:4) {
  spect_man_mo= Mclust(manh_spect_w7[[i]],G = 1:10)
  spect_man_mod_w7 = append(spect_man_mod_w7,spect_man_mo)
  }

spect_man_mod_w9 = numeric(0)
for (i in 1:4) {
  spect_man_mo= Mclust(manh_spect_w9[[i]],G = 1:10)
  spect_man_mod_w9 = append(spect_man_mod_w9,spect_man_mo)
  }

#Mahalanobis

spect_mahal_mod_w7 = numeric(0)
for (i in 1:4) {
  spect_mahal_mo= Mclust(mahal_spect_w7[[i]],G = 1:10)
  spect_mahal_mod_w7 = append(spect_mahal_mod_w7,spect_mahal_mo)
  }

spect_mahal_mod_w9 = numeric(0)
for (i in 1:4) {
  spect_mahal_mo= Mclust(mahal_spect_w9[[i]],G = 1:10)
  spect_mahal_mod_w9 = append(spect_mahal_mod_w9,spect_mahal_mo)
}

# Clustering results
spect_mod_man_w7 = list(spect_man_mod_w7[[14]],spect_man_mod_w7[[29]],spect_man_mod_w7[[44]],spect_man_mod_w7[[59]])
                    
spect_mod_man_w9 = list(spect_man_mod_w9[[14]],spect_man_mod_w9[[29]],spect_man_mod_w9[[44]],spect_man_mod_w9[[59]])
                
spect_mod_mahal_w7 = list(spect_mahal_mod_w7[[14]],spect_mahal_mod_w7[[29]],spect_mahal_mod_w7[[44]],spect_mahal_mod_w7[[59]])                         
spect_mod_mahal_w9 = list(spect_mahal_mod_w9[[14]],spect_mahal_mod_w9[[29]],spect_mahal_mod_w9[[44]],spect_mahal_mod_w9[[59]])
```

Clustering Validation

1. V- Measure validation function - comparison of chemical and spectral classifications

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

2. V-Measure calculation

#Manhattan

# V-Measure Manhattan window size 7

v.measure_manh_w7 = numeric(0)
for (i in 1:4) {
  v = v.measure(chem_mod_man[[i]],spect_mod_man_w7[[i]]) 
  v.measure_manh_w7 = append(v.measure_manh_w7,list(unique(v)))
}

# V-Measure Manhattan window size 9

v.measure_manh_w9 = numeric(0)
for (i in 1:4) {
  v = v.measure(chem_mod_man[[i]],spect_mod_man_w9[[i]]) 
  v.measure_manh_w9 = append(v.measure_manh_w9,list(unique(v)))
}

#Mahalanobis

# V-Measure Mahalanobis window size 7

v.measure_mahal_w7 = numeric(0)
for (i in 1:4) {
  v = v.measure(chem_mod_mahal[[i]],spect_mod_mahal_w7[[i]]) 
  v.measure_mahal_w7 = append(v.measure_mahal_w7,list(unique(v)))
}

# V-Measure Mahalanobis window size 9

v.measure_mahal_w9 = numeric(0)
for (i in 1:4) {
  v = v.measure(chem_mod_mahal[[i]],spect_mod_mahal_w9[[i]]) 
  v.measure_mahal_w9 = append(v.measure_mahal_w9,list(unique(v)))
}

3. V-Measure graph

# Merge v-measure score

v.meas=data.frame(unlist(list(v.measure_manh_w7,v.measure_manh_w9,v.measure_mahal_w7,v.measure_mahal_w9)))

# Create labels

e = c("Manhattan7_w7","Manhattan8_w7","Manhattan9_w7","Manhattan10_w7","Manhattan7_w9","Manhattan8_w9","Manhattan9_w9","Manhattan10_w9",
      "Mahalanobis7_w7","Mahalanobis8_w7","Mahalanobis9_w7","Mahalanobis10_w7","Mahalanobis7_w9","Mahalanobis8_w9","Mahalanobis9_w9","Mahalanobis10_w9")

# combine dataset

v.measure_score = cbind(v.meas,e)
colnames(v.measure_score) = c("V.Measure","Distance_and_window_size")
v.measure_score[,1] = round(v.measure_score[,1],digits =  3)

# Plot V-Measure
ggplot(v.measure_score, aes(x = V.Measure,y = Distance_and_window_size))+
  geom_point(size = 4, col = "blue", shape = 20)+
  scale_x_continuous(breaks=c(seq(0,1,0.1)))+
  geom_text(aes(label=V.Measure,hjust = 0,vjust = 1), size=4, colour="red")+
  ggtitle("V-Measure for Manhattan and Mahalanobis Distance")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size=rel(1.5), face = "bold"))+
  theme(axis.text.x = element_text(face="bold", size=11),
        axis.text.y = element_text(face="bold", size=11))



