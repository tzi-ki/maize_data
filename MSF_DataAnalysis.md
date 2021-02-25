# **Data from B73 and Mo17 maize strains**

## Leaves extracts

[Download the data](<https://doi.org/10.5281/zenodo.4056710>) into your working directory 


```{r}
# Set your working directory
setwd()
```

### Install the libraries needed

```{r eval = FALSE}
if (!require("MALDIquant")) install.packages("MALDIquant")
if (!require("MALDIquantForeign")) install.packages("MALDIquantForeign")
if (!require("pheatmap")) install.packages("pheatmap")
if (!require("factoextra")) install.packages("factoextra")
if (!require("tidyselect")) install.packages("tidyselect")

```

### Load the libraries
```{r message = FALSE}
library(MALDIquant) 
library(MALDIquantForeign) 
library(pheatmap)
library(factoextra)
```

# Mass Spectra Analysis

### View the available files list

```{r}
(mzML_list<- list.files(pattern = "*.mzML"))
```

### Import Data & Sum Scans

```{r warning=FALSE, results='hide'}
filecounter<-1
spectrum<-list()
mzMLno<-length(mzML_list)

while (filecounter<=mzMLno){
  mzML_filename<-mzML_list[filecounter]
  print(mzML_filename)
  importspectrum <- import(mzML_filename)
  sumspectrum <- averageMassSpectra(importspectrum, method="sum")
  spectrum[[filecounter]]<-sumspectrum
  print(filecounter)
  filecounter<-filecounter+1
}
```

### Calibrate intensity using TIC "Total Ion Current" as reference 

```{r}
spectrum_cI_tic <- calibrateIntensity(spectrum, method="TIC")
```

### Plot a mass spectrum 
```{r fig.align='center'}
plot(spectrum_cI_tic[[1]])# Use double square bracket to select an individual spectrum 

```

### Peak picking, binning and peak filtering

```{r results='hide'}
peaks <- detectPeaks(spectrum_cI_tic, SNR = 3)# Adjust signal-to-noise relation 

warpingFunctions <- determineWarpingFunctions(peaks)
peaks <- warpMassPeaks(peaks,warpingFunctions)

peaks <- binPeaks(peaks)
peaks <- filterPeaks(peaks, minFrequency = 0.25)

```

### Generate the Feature Matrix

```{r}
samplenames <- gsub(".mzML","",mzML_list)
featureMatrix <- as.data.frame(intensityMatrix(peaks))

rownames(featureMatrix) <- samplenames
colnames(featureMatrix) <- round(as.numeric(colnames(featureMatrix)),2)

featureMatrix[is.na(featureMatrix)] <- 0
```

If require export the data matrix use `write.csv(featureMatrix,"FeatureMatrix_TIC.csv")` 

# Data visualization

### Find maximal value for each column (mz-bin)

```{r}
columnmax <- apply(featureMatrix, 2, max)
```

### Sort and select the most intense values
```{r}
colmax <- sort(columnmax,decreasing = TRUE)
colmaxselection <- as.matrix(colmax[1:20]) # Use a range to select your ions
colmaxselection <- as.vector(rownames(colmaxselection))

max_normout_featureMatrix <- featureMatrix[,colmaxselection]
t_max_normout_featureMatrix <- t(max_normout_featureMatrix)
```

## Heatmap 

```{r fig.align='center'}
pheatmap(t_max_normout_featureMatrix)# Default graphic: Heatmap arguments,  (distance measure: correlation, and linkage method: complete) 
```

### Personalized Heatmap 
```{r fig.height = 8, fig.align='center'}
pheatmap((t_max_normout_featureMatrix), scale = "row",
         cluster_cols = TRUE, cluster_rows = TRUE,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_methode = "complete", cellwidth = 16, 
         cellheight = 18, fontsize = 10)
```

## Principal Component Analisys (PCA)

### Calculate k-means and visualize cluster assignation
```{r}
km.out <- kmeans(featureMatrix, centers = 2, nstart = 20, iter.max = 50)
km.out$cluster

```

###  Create the PCA
```{r}
pca1 <- prcomp(featureMatrix)
```


### Plot PCA
```{r warning=FALSE, fig.align='center'}
fviz_pca_ind(pca1, axes = c(1,2), 
             col.ind = factor(km.out$cluster),
             addEllipses= T, repel = T,
             legend.title = "Genotype")
```

# Data mining   

## Random Forest

### Create dataframe for the Random Forest model
```{r}
IonFrame <- as.data.frame(featureMatrix)

colnames(IonFrame) <- paste("mz", colnames(IonFrame), sep="")

IonFrame$InbredLine <- as.factor(substr(colnames(t(featureMatrix)), start = 1, stop = 4))# Create variable of interest
```

### *First option*: The Rattle package option

This R package is a Graphical User Interface for Data Science, that enables you to do multiple types of analyzes, Random Forest is one of them.

### Install and load rattle package

```{r eval=FALSE}
if (!require("rattle")) install.packages("rattle")
library(rattle)
```

### Open rattle
```{r eval=FALSE}
rattle() 
```

After the GUI opens go to the `Data tab`, select on the source section the R Dataset option, then choose the option *IonFrame* and click on the execute tab to load it. After loading your data, go to the `Model tab` and select the option *Forest* and click on `Execute tab` to obtain the results.

### *Second option*: The direct option

Without rattle GUI. Using directly the R package `randomForest`

```{r eval = FALSE}
if (!require("randomForest")) install.packages("randomForest")
```
### Generate Random Forest model

```{r message=FALSE}
library(randomForest)
```


```{r}
set.seed(42)
Ionmatrix.rf  <- randomForest(InbredLine ~ ., data=IonFrame, 
                              importance = TRUE, proximity = TRUE)
print(Ionmatrix.rf)

IMP <- as.data.frame(importance(Ionmatrix.rf))# Obtain data frame
IMP <- IMP[order(IMP$MeanDecreaseGini, decreasing = TRUE),]# Order according MeanDecreaseGini
IMP[1:20,]# List the 20 most discriminant ions
```


### Create heatmaps using variables with the highest Gini 
```{r fig.align='center', fig.height=7, fig.width=10}
row.names(IMP) <- substr((row.names(IMP)),3,18) # eliminate mz 
SelIonRF <- rownames(t(featureMatrix)) %in% rownames(IMP[1:20,])# 20 variables

pheatmap(t(featureMatrix[,SelIonRF]), scale = "row",
         cluster_cols = TRUE, cluster_rows = F,
         clustering_distance_cols = "correlation",
         clustering_methode = "complete",cellwidth = 16,
         cellheight = 18, fontsize=10)

```

