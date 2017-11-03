DOG MICROBIOME BASIC ANALYSES
================================
Background
----------

There are three dogs which are treated with increased percentage of a
compound in their diet: 5 different treatments (0-4, representing an
increased percentage of a compound in their diet) Analyses included
here:

-   Import .biom and .tre files generated in QIIME as well as metadata
    file: merge these objects in a phyloseq object
-   Basic data filter: assess number of reads/sample, plot rarefaction
    curves, and exclude low abundance OTUs and samples that do not meet
    minimum number of reads cutoff.
-   Basic exploratory plots including bar plots, alpha- and
    beta-diversity, heatmaps.
-   Differential abundance testing by a) Dog and b) Treatment

## This tutorial should be run in an interactive session. Please do not runnning anything on the headnode.

### To get to a compute node do
```bash
qsub -I -q UCTlong -l walltime=04:00:00
```
Once you are on a compute node you will see that the prompt changes from ```@srvslshpc001``` to ```@srvslshpc60X``` e.g.

```bash
gerrit@srvslshpc001:~> qsub -I -q UCTlong -l walltime=08:00:00
qsub: waiting for job 1598565.srvslshpc001 to start
qsub: job 1598565.srvslshpc001 ready

gerrit@srvslshpc601:~> hostname
srvslshpc601

````

Load the R environment required
-------------------------------

    module load software/R-3.3.0

Startup R 
----------------------

    enter 'R' into terminal


Import data and create phyloseq object
-------------------------------------------

**Import BIOM file (generated in QIIME) into a phyloseq object**


```r
library(phyloseq)
library(ggplot2)
```

```
## Warning: package 'ggplot2' was built under R version 3.3.2
```

```r
library(gridExtra)
library(dunn.test)
library(vegan)
```

```
## Loading required package: permute
```

```
## Loading required package: lattice
```

```
## This is vegan 2.4-0
```

```r
library(randomForest)
```

```
## randomForest 4.6-12
```

```
## Type rfNews() to see new features/changes/bug fixes.
```

```
## 
## Attaching package: 'randomForest'
```

```
## The following object is masked from 'package:gridExtra':
## 
##     combine
```

```
## The following object is masked from 'package:ggplot2':
## 
##     margin
```

```r
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following object is masked from 'package:randomForest':
## 
##     combine
```

```
## The following object is masked from 'package:gridExtra':
## 
##     combine
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```
**Import custom functions used in script**


```r
source("/scratch/DB/bio/training/16SrRNA/16SrRNA-hex-tutorial/src/microbiome_custom_functions.R")
```

```
## Loading required package: pkgmaker
```

```
## Loading required package: registry
```

```
## Loading required package: rngtools
```

```
## Loading required package: cluster
```

```
## NMF - BioConductor layer [OK] | Shared memory capabilities [NO: bigmemory] | Cores 7/8
```

```
##   To enable shared memory capabilities, try: install.extras('
## NMF
## ')
```

```
## 
## Attaching package: 'psych'
```

```
## The following object is masked from 'package:randomForest':
## 
##     outlier
```

```
## The following objects are masked from 'package:ggplot2':
## 
##     %+%, alpha
```

```
## matrixStats v0.50.2 (2016-04-24) successfully loaded. See ?matrixStats for help.
```

```
## 
## Attaching package: 'matrixStats'
```

```
## The following objects are masked from 'package:Biobase':
## 
##     anyMissing, rowMedians
```

```
## The following object is masked from 'package:dplyr':
## 
##     count
```

```
## Loading required package: xtable
```

```
## Loading required package: MASS
```

```
## 
## Attaching package: 'MASS'
```

```
## The following object is masked from 'package:dplyr':
## 
##     select
```

```
## 
## Attaching package: 'fifer'
```

```
## The following object is masked from 'package:Biobase':
## 
##     contents
```

```
## Loading required package: limma
```

```
## 
## Attaching package: 'limma'
```

```
## The following object is masked from 'package:BiocGenerics':
## 
##     plotMA
```

```
## Loading required package: glmnet
```

```
## Loading required package: Matrix
```

```
## Loading required package: foreach
```

```
## foreach: simple, scalable parallel programming from Revolution Analytics
## Use Revolution R for scalability, fault tolerance and more.
## http://www.revolutionanalytics.com
```

```
## Loaded glmnet 2.0-5
```

```
## Loading required package: RColorBrewer
```

```
## Loading required package: gplots
```

```
## 
## Attaching package: 'gplots'
```

```
## The following object is masked from 'package:stats':
## 
##     lowess
```
**Set the working directory and import data**


```r
setwd("/scratch/DB/bio/training/16SrRNA/16SrRNA-hex-tutorial/")
inDir <- getwd()#specify input directory
outDir <- paste0(getwd(),"/results/downstream_analyses") #specify output directory
phy <- import_biom(BIOMfilename = paste0(inDir,"/results/otus_table.tax.biom"), 
		verbose = TRUE)#
ntaxa(phy) #(number of OTUs)
```

```
## [1] 179
```

```r
sample_names(phy) <- sub("\\/1","",sample_names(phy))#remove "/1" from filenames
#add phylogenetic tree (.tre file generated in QIIME)
tree <- read_tree_greengenes(paste0(inDir,"/results/otus_repsetOUT_aligned_pfiltered.tre"))
#merge phy and tree
phy <- merge_phyloseq(phy,tree)
```
**Data cleanup**


```r
colnames(tax_table(phy))
```

```
## [1] "Rank1" "Rank2" "Rank3" "Rank4" "Rank5" "Rank6" "Rank7"
```

```r
colnames(tax_table(phy)) <-  c("Kingdom", "Phylum" , "Class" , "Order" , "Family" , "Genus", "Species")#e.g. replace "Rank1" with "Kingdom"
#clean taxonomic annotations, at the moment they are for example 'k__Bacteria'; 'p_Firmicutes' - remove k__ and p__ ...
tax_table(phy)[,"Kingdom"] <- sub("k__","",tax_table(phy)[,"Kingdom"])
tax_table(phy)[,"Phylum"] <- sub("p__","",tax_table(phy)[,"Phylum"])
tax_table(phy)[,"Class"] <- sub("c__","",tax_table(phy)[,"Class"])
tax_table(phy)[,"Order"] <- sub("o__","",tax_table(phy)[,"Order"])
tax_table(phy)[,"Family"] <- sub("f__","",tax_table(phy)[,"Family"])
tax_table(phy)[,"Genus"] <- sub("g__","",tax_table(phy)[,"Genus"])
tax_table(phy)[,"Species"] <- sub("s__","",tax_table(phy)[,"Species"])
t= which(is.na(tax_table(phy)[,"Phylum"])) 
tax_table(phy) = tax_table(phy)[-t,] #remove rows that don't at least have Phylum-level annotation
```
**Import metadata and merge with phyloseq object**


```r
meta <-  read.table(paste0(inDir,"/practice.dataset1.metadata.tsv"), sep = "\t", header =TRUE, row.names=1)
head(meta)
```

```
##       Dog Treatment
## Dog1    B         2
## Dog2    G         3
## Dog3    K         3
## Dog8    B         4
## Dog9    G         0
## Dog10   K         4
```

```r
rownames(meta)
```

```
##  [1] "Dog1"  "Dog2"  "Dog3"  "Dog8"  "Dog9"  "Dog10" "Dog15" "Dog16"
##  [9] "Dog17" "Dog22" "Dog23" "Dog24" "Dog29" "Dog30" "Dog31"
```

```r
head(sample_names(phy))
```

```
## [1] "Dog10" "Dog15" "Dog16" "Dog17" "Dog1"  "Dog22"
```

```r
length(sample_names(phy))#15
```

```
## [1] 15
```

```r
length(rownames(meta))#15 (check if same number of samples in .biom file and metadatafile)
```

```
## [1] 15
```

```r
length(intersect(rownames(meta),sample_names(phy)))#15 (check that the sample names match in all cases)
```

```
## [1] 15
```

```r
sample_data(phy) <- meta#assign the metadata to the phyloseq object 'phy' (phyloseq will put these in the right order)
nsamples(phy)
```

```
## [1] 15
```

```r
str(sample_data(phy))#need to change treatment column to factor variable
```

```
## 'data.frame':	15 obs. of  2 variables:
## Formal class 'sample_data' [package "phyloseq"] with 4 slots
##   ..@ .Data    :List of 2
##   .. ..$ : Factor w/ 3 levels "B","G","K": 3 1 2 3 1 1 2 3 1 2 ...
##   .. ..$ : int  4 1 4 0 2 3 1 2 0 3 ...
##   ..@ names    : chr  "Dog" "Treatment"
##   ..@ row.names: chr  "Dog10" "Dog15" "Dog16" "Dog17" ...
##   ..@ .S3Class : chr "data.frame"
```

```r
sample_data(phy)[,"Treatment"] <- as.numeric(unlist(sample_data(phy)[,"Treatment"]))
```
**Save phyloseq object as an .RData file**


```r
save(phy, file = paste0(outDir,"/CBIO_16s_cert.RData")) #Save annotated object as a .RData object for quick reload if required at a later stage
#load(paste0(outDir,"/CBIO_16s_cert.RData")) #this is how you would reload the .RData object 'phy'
```
Explore number of reads per sample, make rarefaction curves and filter data as necessary
-------------------------------------------
**Explore number of reads per sample**

```r
reads <- sample_sums(phy)
length(which(reads<5000))
```

```
## [1] 0
```

```r
raremax <- min(reads)
raremax
```

```
## [1] 63980
```

```r
rarecurve(t(otu_table(phy)), step = 100, sample = raremax,xlab = "number of reads/sample", ylab = "number of OTUs",
		label = FALSE, xlim = c(0,100000))
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)

All samples have sufficient sequencing depth for inclusion in downstream analyses. The vertical line in the above plot indicates the sample with the lowest number of reads. Now we will scale data to account for differences in the number of reads/sample and filter rare OTUs that are not of biological interest for the purpose of this analysis (e.g. occurs only in one sample).

**Standardize abundances to median sequence depth**

```r
total = median(sample_sums(phy))
standf = function(x, t=total) round(t * (x / sum(x)))
M.std = transform_sample_counts(phy, standf)
```
**Apply mild OTU filter**

Select OTUs where the rowsum for that OTU has at least 20% of samples with a count of 10 each OR where that OTU > 0.001% of the total median count (for cases where the minority of samples may have high counts of a rare OTU)

```r
M.f = filter_taxa(M.std,function(x) sum(x > 10) > (0.02*length(x)) | sum(x) > 0.001*total, TRUE)
ntaxa(M.f)
```

```
## [1] 135
```
**Basic exploratory plots: alpha- and beta-diversity, barplots, heatmap**
-------------------------------------------
**Alpha diversity by dog**


```r
p <- plot_richness(M.std,x = "Dog",color = "Treatment",measures=c("Shannon"), 
		title = paste0("Standardized to total reads, N=",nsamples(M.std)))+theme(axis.text=element_text(size=16, face="bold"),
				axis.title=element_text(size=16,face="bold"))+geom_point(size=5)
p
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png)

```r
pdf(paste0(outDir,"/alpha_diversity_by_dog_treatment.pdf"))
p
dev.off()
```

```
## quartz_off_screen 
##                 2
```
Is there a significant difference in alpha diversity between dogs irrespective of treatment?

```r
est <- estimate_richness(M.f, split = TRUE, measures = c("Shannon"))
temp <- cbind(est,sample_data(M.f)[,"Dog"])
head(temp)
```

```
##        Shannon Dog
## Dog10 2.859687   K
## Dog15 2.550581   B
## Dog16 2.819618   G
## Dog17 2.823232   K
## Dog1  3.177816   B
## Dog22 2.525255   B
```

```r
t <- kruskal.test(temp[,1]~temp[,2])
t
```

```
## 
## 	Kruskal-Wallis rank sum test
## 
## data:  temp[, 1] by temp[, 2]
## Kruskal-Wallis chi-squared = 2.54, df = 2, p-value = 0.2808
```

```r
dunn.test(temp[,1],temp[,2])#post-hoc testing to see which dogs are different
```

```
##   Kruskal-Wallis rank sum test
## 
## data: x and group
## Kruskal-Wallis chi-squared = 2.54, df = 2, p-value = 0.28
## 
## 
##                            Comparison of x by group                            
##                                 (No adjustment)                                
## Col Mean-|
## Row Mean |          B          G
## ---------+----------------------
##        G |  -1.343502
##          |     0.0896
##          |
##        K |   0.070710   1.414213
##          |     0.4718     0.0786
```
Dog G has higher alpha diversity than dogs K and B irrespective of treatment, but this difference is not significant

**Alpha diversity by treatment**

```r
p <- plot_richness(M.std,x = "Treatment",color = "Dog",measures=c("Shannon"), 
				title = paste0("Standardized to total reads, N=",nsamples(M.std)))+theme(axis.text=element_text(size=16, face="bold"),
				axis.title=element_text(size=16,face="bold"))+geom_point(size=5)
p
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-1.png)

```r
pdf(paste0(outDir,"/alpha_diversity_by_treatment_dog.pdf"))
p
dev.off()
```

```
## quartz_off_screen 
##                 2
```
Are there significant differences in alpha diversity by treatment?

```r
temp <- cbind(est,sample_data(M.f)[,"Treatment"])
head(temp)
```

```
##        Shannon Treatment
## Dog10 2.859687         4
## Dog15 2.550581         1
## Dog16 2.819618         4
## Dog17 2.823232         0
## Dog1  3.177816         2
## Dog22 2.525255         3
```

```r
t <- kruskal.test(temp[,1]~temp[,2])
t
```

```
## 
## 	Kruskal-Wallis rank sum test
## 
## data:  temp[, 1] by temp[, 2]
## Kruskal-Wallis chi-squared = 1.5667, df = 4, p-value = 0.8148
```

```r
dunn.test(temp[,1],temp[,2])
```

```
##   Kruskal-Wallis rank sum test
## 
## data: x and group
## Kruskal-Wallis chi-squared = 1.5667, df = 4, p-value = 0.81
## 
## 
##                            Comparison of x by group                            
##                                 (No adjustment)                                
## Col Mean-|
## Row Mean |          4          1          0          2
## ---------+--------------------------------------------
##        1 |   0.365148
##          |     0.3575
##          |
##        0 |  -0.365148  -0.730296
##          |     0.3575     0.2326
##          |
##        2 |  -0.091287  -0.456435   0.273861
##          |     0.4636     0.3240     0.3921
##          |
##        3 |  -0.821583  -1.186732  -0.456435  -0.730296
##          |     0.2057     0.1177     0.3240     0.2326
```

**Beta diversity using NMDS with Bray-Curtis as distance metric**

```r
set.seed(2)
GP.ord.BC <- ordinate(M.f, "NMDS", "bray", k=2, trymax=100)
```

```
## Square root transformation
## Wisconsin double standardization
## Run 0 stress 0.07152504 
## Run 1 stress 0.07375128 
## Run 2 stress 0.074195 
## Run 3 stress 0.1164627 
## Run 4 stress 0.07282976 
## Run 5 stress 0.3620394 
## Run 6 stress 0.07419498 
## Run 7 stress 0.07327274 
## Run 8 stress 0.2139577 
## Run 9 stress 0.07375192 
## Run 10 stress 0.07152511 
## ... Procrustes: rmse 3.59076e-05  max resid 0.0001198476 
## ... Similar to previous best
## Run 11 stress 0.07348428 
## Run 12 stress 0.07282963 
## Run 13 stress 0.07419507 
## Run 14 stress 0.0732727 
## Run 15 stress 0.07152527 
## ... Procrustes: rmse 0.0001508438  max resid 0.0003573785 
## ... Similar to previous best
## Run 16 stress 0.1160219 
## Run 17 stress 0.07375338 
## Run 18 stress 0.07419527 
## Run 19 stress 0.07152524 
## ... Procrustes: rmse 5.606673e-05  max resid 0.0001361385 
## ... Similar to previous best
## Run 20 stress 0.1164624 
## *** Solution reached
```

```r
GP.ord.BC
```

```
## 
## Call:
## metaMDS(comm = veganifyOTU(physeq), distance = distance, k = 2,      trymax = 100) 
## 
## global Multidimensional Scaling using monoMDS
## 
## Data:     wisconsin(sqrt(veganifyOTU(physeq))) 
## Distance: bray 
## 
## Dimensions: 2 
## Stress:     0.07152504 
## Stress type 1, weak ties
## Two convergent solutions found after 20 tries
## Scaling: centring, PC rotation, halfchange scaling 
## Species: expanded scores based on 'wisconsin(sqrt(veganifyOTU(physeq)))'
```

```r
color = c("Treatment")
shape = c("Dog")
title=c("NMDS of 16S microbiome,Bray-Curtis distance,k=2")
MDS = plot_ordination(M.f, GP.ord.BC, color = color,shape=shape, 
		title = title)
MDS.1  = MDS +theme(axis.text=element_text(size=16, face="bold"),
				axis.title=element_text(size=18,face="bold"), legend.title=element_text(size=14))+
		theme_bw()+labs(color=color, shape=shape)+geom_point(size=5)

MDS.1
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-1.png)

```r
pdf(paste0(outDir,"/NMDS_Dogs_tretment_Bray_Curtis.pdf"),8,5)
MDS.1
dev.off()
```

```
## quartz_off_screen 
##                 2
```
**Beta diversity using NMDS with Unifrac as distance metric**

```r
GP.ord.U <- ordinate(M.f, "NMDS", "unifrac")
```

```
## Run 0 stress 0.09882137 
## Run 1 stress 0.1061522 
## Run 2 stress 0.09882136 
## ... New best solution
## ... Procrustes: rmse 4.127018e-05  max resid 8.351313e-05 
## ... Similar to previous best
## Run 3 stress 0.09882146 
## ... Procrustes: rmse 4.510627e-05  max resid 9.586255e-05 
## ... Similar to previous best
## Run 4 stress 0.3349044 
## Run 5 stress 0.1100681 
## Run 6 stress 0.09925322 
## ... Procrustes: rmse 0.05419579  max resid 0.1641213 
## Run 7 stress 0.1100682 
## Run 8 stress 0.1100683 
## Run 9 stress 0.2113797 
## Run 10 stress 0.2258924 
## Run 11 stress 0.09882136 
## ... Procrustes: rmse 1.99549e-05  max resid 3.926798e-05 
## ... Similar to previous best
## Run 12 stress 0.1348672 
## Run 13 stress 0.1672158 
## Run 14 stress 0.1100681 
## Run 15 stress 0.2439499 
## Run 16 stress 0.1100685 
## Run 17 stress 0.1061522 
## Run 18 stress 0.1583075 
## Run 19 stress 0.09925322 
## ... Procrustes: rmse 0.0541974  max resid 0.164135 
## Run 20 stress 0.09882146 
## ... Procrustes: rmse 0.0001361423  max resid 0.0002761411 
## ... Similar to previous best
## *** Solution reached
```

```r
GP.ord.U
```

```
## 
## Call:
## metaMDS(comm = ps.dist) 
## 
## global Multidimensional Scaling using monoMDS
## 
## Data:     ps.dist 
## Distance: user supplied 
## 
## Dimensions: 2 
## Stress:     0.09882136 
## Stress type 1, weak ties
## Two convergent solutions found after 20 tries
## Scaling: centring, PC rotation 
## Species: scores missing
```

```r
color = c("Treatment")
shape = c("Dog")

title=c("NMDS of 16S microbiome, Unifrac distance, k=2")

MDS = plot_ordination(M.f, GP.ord.U, color = color, shape=shape, 
		title = title)
MDS.1  = MDS +theme(axis.text=element_text(size=16, face="bold"),
				axis.title=element_text(size=18,face="bold"), legend.title=element_text(size=14))+
		theme_bw()+labs(color=color)+geom_point(size=5)
MDS.1
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15-1.png)

```r
pdf(paste0(outDir,"/NMDS_Dogs_tretment_Bray_Curtis.pdf"),8,5)
MDS.1
dev.off()
```

```
## quartz_off_screen 
##                 2
```
**Create a heatmap of taxa merged at the lowest available taxonomic level**

```r
M.phy <- tax_glom.kv(M.f)#this function is available in the 'microbiome_custom_functions.R' script loaded at the beginning of this script
```

```
## [1] "Removing phylogenetic tree"
## [1] "There are now 58 merged taxa"
```

```r
ntaxa(M.phy)
```

```
## [1] 58
```

```r
filename <- c("cbio_cert_heatmap_merged_taxa")
main <- paste("Merged taxa, Bray-Curtis distance")
f = paste0(outDir,"/",filename,".pdf")
#color specification for column annotations above heatmap:
D.cols = c("B"="#CC79A7","G"="#56B4E9","K"="#F0E442")
colours = list(Dog=D.cols)

#create distance matrix and calculate tree:
set.seed(2)
diss <- distance(M.phy,method = "bray", type = "samples")
clust.res<-hclust(diss)
sample.order = clust.res$order
#heatmap is output to file (the heatmap.k function can be found in the 'microbiome_custom_functions.R' script)
hm = heatmap.k(physeq= M.phy,
		annot.cols = c(1,2),
		main = main,filename = f,colours=colours,Colv = sample.order,labrow = TRUE)	
```

```
## [1] "including all samples"
## [1] "including all otus"
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16-1.png)

```r
print(hm)
```

```
## $Rowv
## 'dendrogram' with 2 branches and 58 members total, at height 31.65134 
## 
## $rowInd
##  [1] 54 22 19 28 44 30 33 35 17 37 39 56 51 57 20  4 50 42 25  8 40  5 21
## [24] 15  3 49 13 10 18 43 27 26 14  7  6  9 52 55 16 34 41 47 38 12 48 11
## [47] 53 45  2 36 24 31 58 23 29 32 46  1
## 
## $colInd
##  [1] 13  9  2  6  5 10 15  3  7 11  1 14 12  4  8
## 
## $col
##     #4575B4     #4979B6     #4E7DB8     #5282BB     #5786BD     #5C8BBF 
## -3.34029373 -3.15627003 -2.97224633 -2.78822263 -2.60419892 -2.42017522 
##     #608FC2     #6594C4     #6998C6     #6E9DC9     #73A1CB     #77A6CE 
## -2.23615152 -2.05212782 -1.86810411 -1.68408041 -1.50005671 -1.31603301 
##     #7CAAD0     #80AFD2     #85B3D5     #8AB8D7     #8EBCD9     #93C0DB 
## -1.13200931 -0.94798560 -0.76396190 -0.57993820 -0.39591450 -0.21189079 
##     #98C3DD     #9CC6DF     #A1CAE1     #A6CDE2     #ABD0E4     #B0D3E6 
## -0.02786709  0.15615661  0.34018031  0.52420402  0.70822772  0.89225142 
##     #B4D6E8     #B9D9E9     #BEDCEB     #C3E0ED     #C8E3EF     #CCE6F0 
##  1.07627512  1.26029883  1.44432253  1.62834623  1.81236993  1.99639363 
##     #D1E9F2     #D6ECF4     #DBEFF6     #E0F3F7     #E1F3F4     #E3F4F1 
##  2.18041734  2.36444104  2.54846474  2.73248844  2.91651215  3.10053585 
##     #E5F5ED     #E7F5EA     #E9F6E6     #EBF7E3     #EDF8DF     #EFF8DC 
##  3.28455955  3.46858325  3.65260696  3.83663066  4.02065436  4.20467806 
##     #F0F9D8     #F2FAD5     #F4FBD2     #F6FBCE     #F8FCCB     #FAFDC7 
##  4.38870177  4.57272547  4.75674917  4.94077287  5.12479657  5.30882028 
##     #FCFDC4     #FEFEC0     #FEFEBD     #FEFCBA     #FEFAB7     #FEF8B5 
##  5.49284398  5.67686768  5.86089138  6.04491509  6.22893879  6.41296249 
##     #FEF6B2     #FEF4AF     #FEF2AC     #FEF0A9     #FEEFA6     #FEEDA3 
##  6.59698619  6.78100990  6.96503360  7.14905730  7.33308100  7.51710471 
##     #FEEBA1     #FEE99E     #FEE79B     #FEE598     #FEE395     #FEE192 
##  7.70112841  7.88515211  8.06917581  8.25319951  8.43722322  8.62124692 
##     #FEDF8F     #FDDA8C     #FDD589     #FDD085     #FDCB82     #FDC67F 
##  8.80527062  8.98929432  9.17331803  9.35734173  9.54136543  9.72538913 
##     #FDC17B     #FDBC78     #FDB775     #FCB271     #FCAD6E     #FCA86B 
##  9.90941284 10.09343654 10.27746024 10.46148394 10.64550765 10.82953135 
##     #FCA367     #FC9E64     #FC9961     #FC945D     #FC8F5A     #FA8A57 
## 11.01355505 11.19757875 11.38160246 11.56562616 11.74964986 11.93367356 
##     #F88454     #F67E51     #F4794E     #F1734B     #EF6D48     #ED6845 
## 12.11769726 12.30172097 12.48574467 12.66976837 12.85379207 13.03781578 
##     #EB6242     #E85D3F     #E6573C     #E45139     #E24C36     #DF4633 
## 13.22183948 13.40586318 13.58988688 13.77391059 13.95793429 14.14195799 
##     #DD4030     #DB3B2D     #D9352A     #D73027        <NA> 
## 14.32598169 14.51000540 14.69402910 14.87805280 15.06207650
```
**Barplots by dog**
------------------------------

```r
level = "Genus"
count = 500
perc = 0.25
#barplot will be written to file (the bar.plots function can be found in the 'microbiome_custom_functions.R' script)
barplot = bar.plots(physeq = M.std,cat = "Dog",level = level, count = count, perc = perc, outDir=outDir, 
		filen = 'Barplots_by_Dog')
print(barplot)
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-17-1.png)

Detect taxa/OTUs that differ significantly by Dog
-------------------------------------------
convert phyloseq object to metagenomeSeq obj. NB use raw data not standardized:

```r
Mraw.f = filter_taxa(phy,function(x) sum(x > 10) > (0.02*length(x)) | sum(x) > 0.001*total, TRUE)
ntaxa(Mraw.f)
```

```
## [1] 135
```

```r
MGS=make_metagenomeSeq(Mraw.f)
```

```
## Default value being used.
```

```r
MGS
```

```
## MRexperiment (storageMode: environment)
## assayData: 135 features, 15 samples 
##   element names: counts 
## protocolData: none
## phenoData
##   sampleNames: Dog10 Dog15 ... Dog9 (15 total)
##   varLabels: Dog Treatment
##   varMetadata: labelDescription
## featureData
##   featureNames: OTU_30 OTU_41 ... OTU_66 (135 total)
##   fvarLabels: OTUname
##   fvarMetadata: labelDescription
## experimentData: use 'experimentData(object)'
## Annotation:
```

**Use Random forests analysis to detect taxa that are good predictors of Dog**

Example used: Dog G vs. Dog B (all treatment points)


```r
sub.index <- sample_names(M.f)[sample_data(M.f)[,"Dog"] != "K"]
phy.temp <- prune_samples(sub.index, M.f)
nsamples(phy.temp)
```

```
## [1] 10
```

```r
RF.k(data = phy.temp, var = "Dog", ntree=10000, cv.fold=10, outDir = outDir, Nfeatures.validation = 3)
```

```
## [1] "0 samples did not have response variable data, removing these..."
## [1] "Data set size:  10 samples with 5 and 5 samples per class"
## [1] "Cross-validated error rates associated with stepwise reduction of features:"
## 135  68  34  17   8   4   2   1 
## 0.0 0.0 0.0 0.0 0.0 0.0 0.2 0.0
```

```
## [1] "*****************************"
## [1] "THE TOP 20 MOST IMPORTANT FEATURES WERE:"
##         predictors         B         G MeanDecreaseAccuracy
## OTU_121    OTU_121  9.671867  9.681629             9.993612
## OTU_51      OTU_51  8.804326  8.794635             9.104466
## OTU_8        OTU_8 10.147888 10.256865            10.526653
## OTU_72      OTU_72  9.560636  9.639463            10.031998
## OTU_104    OTU_104  9.232706  9.125792             9.443961
## OTU_58      OTU_58  9.391857  9.411628             9.868075
## OTU_11      OTU_11  8.359543  8.290286             8.680334
## OTU_119    OTU_119  9.141466  9.134742             9.416359
## OTU_115    OTU_115  8.950253  8.937649             9.235798
## OTU_53      OTU_53  8.925435  8.801074             9.192012
## OTU_74      OTU_74  9.342631  9.406728             9.785011
## OTU_59      OTU_59  8.384654  8.370093             8.661476
## OTU_73      OTU_73  8.405341  8.385556             8.701166
## OTU_140    OTU_140  8.380246  8.371369             8.789388
## OTU_96      OTU_96  7.817527  7.903793             8.080081
## OTU_78      OTU_78  7.751269  7.708754             8.008854
## OTU_150    OTU_150  8.586375  8.567680             8.865167
## OTU_24      OTU_24  7.428018  7.468994             7.720773
## OTU_10      OTU_10  9.535855  9.464636             9.827861
## OTU_29      OTU_29  8.973883  8.878966             9.301188
##         MeanDecreaseGini                 tax
## OTU_121       0.08080000             P.copri
## OTU_51        0.08008667       Clostridiales
## OTU_8         0.07962000           P.copri.1
## OTU_72        0.07860000       Adlercreutzia
## OTU_104       0.07852000         Odoribacter
## OTU_58        0.07734000  Anaerobiospirillum
## OTU_11        0.07688000  Enterobacteriaceae
## OTU_119       0.07652000          Sutterella
## OTU_115       0.07576000        Enterococcus
## OTU_53        0.07570667      Clostridiaceae
## OTU_74        0.07478000  [Mogibacteriaceae]
## OTU_59        0.07470000      [Ruminococcus]
## OTU_73        0.07414000 Erysipelotrichaceae
## OTU_140       0.07386000       Fusobacterium
## OTU_96        0.07374000        Oscillospira
## OTU_78        0.07370000         Bacteroides
## OTU_150       0.07352667     Lachnospiraceae
## OTU_24        0.07346000   Lachnospiraceae.1
## OTU_10        0.07342000          B.producta
## OTU_29        0.07286000        Turicibacter
## [1] "*****************************"
```

```
## [1] "Training AUC=1"
## [1] "Training PPV=1"
## [1] "Training NPV=1"
```

```
## [1] "*****************************"
## [1] "Training set classification summary if using the top 3 features only"
## [1] "Feature(s) selected: OTU_121" "Feature(s) selected: OTU_51" 
## [3] "Feature(s) selected: OTU_8"  
##                 B         G MeanDecreaseAccuracy MeanDecreaseGini
## OTU_8   0.1447083 0.1438417            0.1290467          1.52652
## OTU_121 0.1430083 0.1389750            0.1261224          1.50398
## OTU_51  0.1139500 0.1111500            0.1011867          1.47938
## 
## Call:
##  randomForest(formula = response ~ ., data = rf.data[, c(goodPredictors,      "response")], importance = T, proximity = T, ntree = ntree,      na.action = na.omit) 
##                Type of random forest: classification
##                      Number of trees: 10000
## No. of variables tried at each split: 1
## 
##         OOB estimate of  error rate: 0%
## Confusion matrix:
##   B G class.error
## B 5 0           0
## G 0 5           0
```

The class error rates are 0% (even one OTU enough to discriminate between Dog G and B?)

What if we used merged OTUs?

```r
merged.phy <- tax_glom.kv(phy.temp)
```

```
## [1] "Removing phylogenetic tree"
## [1] "There are now 58 merged taxa"
```

```r
RF.k(data = merged.phy, var = "Dog", ntree=10000, cv.fold=10, outDir = outDir, Nfeatures.validation = 3, descriptor = "merged_OTUs") #for details on RF.k() see microbiome_custom_functions.R file
```

```
## [1] "0 samples did not have response variable data, removing these..."
## [1] "Data set size:  10 samples with 5 and 5 samples per class"
## [1] "Cross-validated error rates associated with stepwise reduction of features:"
##  58  29  14   7   4   1 
## 0.0 0.0 0.0 0.0 0.0 0.1
```

```
## [1] "*****************************"
## [1] "THE TOP 20 MOST IMPORTANT FEATURES WERE:"
##         predictors        B        G MeanDecreaseAccuracy MeanDecreaseGini
## OTU_59      OTU_59 12.84573 12.75514             13.22767        0.1701578
## OTU_27      OTU_27 14.14560 14.15688             14.66689        0.1664900
## OTU_104    OTU_104 13.79309 13.76287             14.34825        0.1652714
## OTU_68      OTU_68 12.54870 12.48715             13.03041        0.1639057
## OTU_72      OTU_72 13.67522 13.74857             14.21517        0.1631848
## OTU_76      OTU_76 13.48093 13.57423             14.07012        0.1623467
## OTU_29      OTU_29 13.77678 13.69971             14.25054        0.1618883
## OTU_115    OTU_115 13.78611 13.66527             14.26590        0.1615750
## OTU_98      OTU_98 13.23068 13.23404             13.60197        0.1611457
## OTU_57      OTU_57 12.79080 12.78306             13.28220        0.1607600
## OTU_26      OTU_26 12.70062 12.54606             13.09699        0.1593200
## OTU_8        OTU_8 13.81179 13.67807             14.29226        0.1565095
## OTU_11      OTU_11 12.19598 12.23061             12.70349        0.1563300
## OTU_55      OTU_55 13.95863 13.71404             14.38892        0.1548617
## OTU_12      OTU_12 11.50679 11.36390             11.94686        0.1536400
## OTU_35      OTU_35 12.44178 12.45397             13.04827        0.1484333
## OTU_74      OTU_74 13.36003 13.41300             13.80032        0.1483867
## OTU_19      OTU_19 13.36975 13.20095             13.78641        0.1470000
## OTU_10      OTU_10 12.36814 12.38032             12.83683        0.1421733
## OTU_61      OTU_61 12.03922 12.09042             12.51234        0.1395000
##                           tax
## OTU_59         [Ruminococcus]
## OTU_27           Ruminococcus
## OTU_104           Odoribacter
## OTU_68     Anaerobiospirillum
## OTU_72          Adlercreutzia
## OTU_76           Oscillospira
## OTU_29           Turicibacter
## OTU_115          Enterococcus
## OTU_98                  Dorea
## OTU_57    Succinivibrionaceae
## OTU_26                  S24-7
## OTU_8                 P.copri
## OTU_11     Enterobacteriaceae
## OTU_55          Bacteroidales
## OTU_12          F.prausnitzii
## OTU_35            Allobaculum
## OTU_74     [Mogibacteriaceae]
## OTU_19               R.gnavus
## OTU_10             B.producta
## OTU_61  Peptostreptococcaceae
## [1] "*****************************"
```

```
## [1] "Training AUC=1"
## [1] "Training PPV=1"
## [1] "Training NPV=1"
```

```
## [1] "*****************************"
## [1] "Training set classification summary if using the top 3 features only"
## [1] "Feature(s) selected: OTU_59"  "Feature(s) selected: OTU_27" 
## [3] "Feature(s) selected: OTU_104"
##                 B         G MeanDecreaseAccuracy MeanDecreaseGini
## OTU_104 0.1392500 0.1415667            0.1262033          1.52932
## OTU_27  0.1363083 0.1369833            0.1223476          1.49262
## OTU_59  0.1169167 0.1174750            0.1051862          1.48760
## 
## Call:
##  randomForest(formula = response ~ ., data = rf.data[, c(goodPredictors,      "response")], importance = T, proximity = T, ntree = ntree,      na.action = na.omit) 
##                Type of random forest: classification
##                      Number of trees: 10000
## No. of variables tried at each split: 1
## 
##         OOB estimate of  error rate: 0%
## Confusion matrix:
##   B G class.error
## B 5 0           0
## G 0 5           0
```

```r
#Note that Nfeatures.validation: 'x' number of top taxa to test (e.g. how good are the top 3 most important taxa at classifying)
```
**Differential abundance testing using MetagenomeSeq package**

Lets again compare dog G vs. dog B (merged taxa), this time using differential abundance testing


```r
colours = list(Dog=D.cols)
a = super.fitZig.kv(physeq = merged.phy,factor = "Dog",outDir = outDir,FileName =c("1_25FC_0.2_Dog_GvsB_taxa_merged"),
		heatmap.descriptor=c("tax_annot"), main=c("Dog G vs. B, taxa merged"), subt=c("subt = FDR < 0.05,|coeff| >= 1.25, >20%+ in either group"), 
		ordered=TRUE, p=0.05, FC = 1.25, perc=0.2, extra.cols = c("Treatment"))
```

```
## [1] "0 of 10 samples were removed due to missing data"
```

```
## Default value being used.
```

```
## [1] "Dog will be modeled as a binary categorical predictor variable"
```

```
## Default value being used.
```

```
## it= 0, nll=18.47, log10(eps+1)=Inf, stillActive=58
## it= 1, nll=19.48, log10(eps+1)=0.01, stillActive=4
## it= 2, nll=19.52, log10(eps+1)=0.03, stillActive=1
## it= 3, nll=19.55, log10(eps+1)=0.01, stillActive=1
## it= 4, nll=19.63, log10(eps+1)=0.00, stillActive=0
## There were  25 OTUs significantly different between B vs. G that met 
##  threshold criteria of p 0.05 absolute FC 1.25 and percentage presence in at least one group of 20 % 
## [1] "writing results and model to file"
## [1] "/scratch/DB/bio/training/16SrRNA/16SrRNA-hex-tutorial/results/downstream_analyses/1_25FC_0.2_Dog_GvsB_taxa_merged_tax_annot.pdf"
```

![plot of chunk unnamed-chunk-21](figure/unnamed-chunk-21-1.png)![plot of chunk unnamed-chunk-21](figure/unnamed-chunk-21-2.png)

```
## [1] "making heatmap of results"
```

```r
print(a)
```

```
##         percent_positive_group0 percent_positive_group1
## OTU_11                      100                     100
## OTU_111                      80                      20
## OTU_19                      100                     100
## OTU_77                      100                      80
## OTU_29                      100                     100
## OTU_115                     100                      40
## OTU_61                      100                     100
## OTU_98                      100                     100
## OTU_10                      100                     100
## OTU_59                      100                     100
## OTU_12                      100                     100
## OTU_73                      100                     100
## OTU_44                       80                     100
## OTU_76                      100                     100
## OTU_27                      100                     100
## OTU_74                       60                     100
## OTU_69                       60                     100
## OTU_39                      100                     100
## OTU_8                       100                     100
## OTU_26                      100                     100
## OTU_57                       60                     100
## OTU_72                        0                     100
## OTU_68                       40                     100
## OTU_55                       20                     100
## OTU_36                       40                     100
##         +samples in group 0 +samples in group 1 mean_positive_group0
## OTU_11                    5                   5                 3645
## OTU_111                   4                   1                   32
## OTU_19                    5                   5                 1142
## OTU_77                    5                   4                   66
## OTU_29                    5                   5                 1147
## OTU_115                   5                   2                   10
## OTU_61                    5                   5                  160
## OTU_98                    5                   5                 1570
## OTU_10                    5                   5                 4032
## OTU_59                    5                   5                  105
## OTU_12                    5                   5                 1821
## OTU_73                    5                   5                  119
## OTU_44                    4                   5                   37
## OTU_76                    5                   5                   16
## OTU_27                    5                   5                  242
## OTU_74                    3                   5                    7
## OTU_69                    3                   5                    4
## OTU_39                    5                   5                   85
## OTU_8                     5                   5                 1524
## OTU_26                    5                   5                   20
## OTU_57                    3                   5                    2
## OTU_72                    0                   5                  NaN
## OTU_68                    2                   5                    3
## OTU_55                    1                   5                    4
## OTU_36                    2                   5                   10
##         mean_positive_group1 oddsRatio      lower       upper     fisherP
## OTU_11                    70    0.0000 0.00000000         Inf 1.000000000
## OTU_111                    1   10.9072 0.45473998 968.7617574 0.206349206
## OTU_19                    85    0.0000 0.00000000         Inf 1.000000000
## OTU_77                     6       Inf 0.02564066         Inf 1.000000000
## OTU_29                   241    0.0000 0.00000000         Inf 1.000000000
## OTU_115                    1       Inf 0.49337123         Inf 0.166666667
## OTU_61                    27    0.0000 0.00000000         Inf 1.000000000
## OTU_98                   409    0.0000 0.00000000         Inf 1.000000000
## OTU_10                  1277    0.0000 0.00000000         Inf 1.000000000
## OTU_59                    32    0.0000 0.00000000         Inf 1.000000000
## OTU_12                   588    0.0000 0.00000000         Inf 1.000000000
## OTU_73                    38    0.0000 0.00000000         Inf 1.000000000
## OTU_44                   201    0.0000 0.00000000  39.0005500 1.000000000
## OTU_76                   121    0.0000 0.00000000         Inf 1.000000000
## OTU_27                  1453    0.0000 0.00000000         Inf 1.000000000
## OTU_74                    73    0.0000 0.00000000   5.1183766 0.444444444
## OTU_69                    88    0.0000 0.00000000   5.1183766 0.444444444
## OTU_39                   598    0.0000 0.00000000         Inf 1.000000000
## OTU_8                  15540    0.0000 0.00000000         Inf 1.000000000
## OTU_26                   456    0.0000 0.00000000         Inf 1.000000000
## OTU_57                   104    0.0000 0.00000000   5.1183766 0.444444444
## OTU_72                    52    0.0000 0.00000000   0.4353226 0.007936508
## OTU_68                   231    0.0000 0.00000000   2.0268713 0.166666667
## OTU_55                   264    0.0000 0.00000000   0.9757790 0.047619048
## OTU_36                   782    0.0000 0.00000000   2.0268713 0.166666667
##         fisherAdjP     coeff      pvalues  adjPvalues  Kingdom
## OTU_11   1.0000000 -5.584733 2.632486e-03 0.008981422 Bacteria
## OTU_111  1.0000000 -3.821533 1.375479e-03 0.004986112 Bacteria
## OTU_19   1.0000000 -3.737026 2.629153e-04 0.001694343 Bacteria
## OTU_77   1.0000000 -3.579532 4.548917e-03 0.013191859 Bacteria
## OTU_29   1.0000000 -2.453196 6.035145e-04 0.002936905 Bacteria
## OTU_115  1.0000000 -2.379397 5.497207e-03 0.015182763 Bacteria
## OTU_61   1.0000000 -2.225829 6.076356e-04 0.002936905 Bacteria
## OTU_98   1.0000000 -2.090437 1.479472e-04 0.001694343 Bacteria
## OTU_10   1.0000000 -1.715693 1.274277e-03 0.004927205 Bacteria
## OTU_59   1.0000000 -1.692829 9.068161e-04 0.004045795 Bacteria
## OTU_12   1.0000000 -1.610710 3.120434e-03 0.010054733 Bacteria
## OTU_73   1.0000000 -1.608000 9.476365e-03 0.023896920 Bacteria
## OTU_44   1.0000000  1.851675 1.996187e-02 0.044530317 Bacteria
## OTU_76   1.0000000  3.097858 5.065905e-04 0.002936905 Bacteria
## OTU_27   1.0000000  3.166680 3.369510e-03 0.010285872 Bacteria
## OTU_74   1.0000000  3.576528 2.139827e-04 0.001694343 Bacteria
## OTU_69   1.0000000  3.815334 1.891429e-02 0.043881157 Bacteria
## OTU_39   1.0000000  3.980148 6.443652e-03 0.016987809 Bacteria
## OTU_8    1.0000000  4.115082 1.219577e-03 0.004927205 Bacteria
## OTU_26   1.0000000  4.973902 2.611733e-04 0.001694343 Bacteria
## OTU_57   1.0000000  5.136267 1.493929e-04 0.001694343 Bacteria
## OTU_72   0.4603175  5.546255 2.124125e-04 0.001694343 Bacteria
## OTU_68   1.0000000  6.085491 1.169281e-04 0.001694343 Bacteria
## OTU_55   0.5523810  6.572744 9.392677e-05 0.001694343 Bacteria
## OTU_36   1.0000000  6.638985 1.950259e-04 0.001694343 Bacteria
##                 Phylum               Class              Order
## OTU_11  Proteobacteria Gammaproteobacteria  Enterobacteriales
## OTU_111     Firmicutes          Clostridia      Clostridiales
## OTU_19      Firmicutes          Clostridia      Clostridiales
## OTU_77      Firmicutes     Erysipelotrichi Erysipelotrichales
## OTU_29      Firmicutes             Bacilli   Turicibacterales
## OTU_115     Firmicutes             Bacilli    Lactobacillales
## OTU_61      Firmicutes          Clostridia      Clostridiales
## OTU_98      Firmicutes          Clostridia      Clostridiales
## OTU_10      Firmicutes          Clostridia      Clostridiales
## OTU_59      Firmicutes          Clostridia      Clostridiales
## OTU_12      Firmicutes          Clostridia      Clostridiales
## OTU_73      Firmicutes     Erysipelotrichi Erysipelotrichales
## OTU_44   Bacteroidetes         Bacteroidia      Bacteroidales
## OTU_76      Firmicutes          Clostridia      Clostridiales
## OTU_27      Firmicutes          Clostridia      Clostridiales
## OTU_74      Firmicutes          Clostridia      Clostridiales
## OTU_69   Bacteroidetes         Bacteroidia      Bacteroidales
## OTU_39   Bacteroidetes         Bacteroidia      Bacteroidales
## OTU_8    Bacteroidetes         Bacteroidia      Bacteroidales
## OTU_26   Bacteroidetes         Bacteroidia      Bacteroidales
## OTU_57  Proteobacteria Gammaproteobacteria      Aeromonadales
## OTU_72  Actinobacteria      Coriobacteriia   Coriobacteriales
## OTU_68  Proteobacteria Gammaproteobacteria      Aeromonadales
## OTU_55   Bacteroidetes         Bacteroidia      Bacteroidales
## OTU_36   Bacteroidetes         Bacteroidia      Bacteroidales
##                        Family              Genus     Species
## OTU_11     Enterobacteriaceae               <NA>        <NA>
## OTU_111       Lachnospiraceae       Epulopiscium        <NA>
## OTU_19        Lachnospiraceae     [Ruminococcus]      gnavus
## OTU_77    Erysipelotrichaceae      [Eubacterium]    dolichum
## OTU_29      Turicibacteraceae       Turicibacter        <NA>
## OTU_115       Enterococcaceae       Enterococcus        <NA>
## OTU_61  Peptostreptococcaceae               <NA>        <NA>
## OTU_98        Lachnospiraceae              Dorea        <NA>
## OTU_10        Lachnospiraceae            Blautia    producta
## OTU_59        Lachnospiraceae     [Ruminococcus]        <NA>
## OTU_12        Ruminococcaceae   Faecalibacterium prausnitzii
## OTU_73    Erysipelotrichaceae               <NA>        <NA>
## OTU_44         Bacteroidaceae        Bacteroides coprophilus
## OTU_76        Ruminococcaceae       Oscillospira        <NA>
## OTU_27        Ruminococcaceae       Ruminococcus        <NA>
## OTU_74     [Mogibacteriaceae]               <NA>        <NA>
## OTU_69         Bacteroidaceae        Bacteroides   uniformis
## OTU_39     Porphyromonadaceae    Parabacteroides        <NA>
## OTU_8          Prevotellaceae         Prevotella       copri
## OTU_26                  S24-7               <NA>        <NA>
## OTU_57    Succinivibrionaceae               <NA>        <NA>
## OTU_72      Coriobacteriaceae      Adlercreutzia        <NA>
## OTU_68    Succinivibrionaceae Anaerobiospirillum        <NA>
## OTU_55                   <NA>               <NA>        <NA>
## OTU_36   [Paraprevotellaceae]               <NA>        <NA>
```
Now again compare dog G vs. dog B with differential abundance testing (individual OTUs as opposed to merged this time)


```r
b = super.fitZig.kv(physeq = phy.temp,factor = "Dog",outDir = outDir,FileName =c("1_25FC_0.2_Dog_GvsB_OTUs"),
		heatmap.descriptor=c("tax_annot"), main=c("Dog G vs. B, OTUs"), subt=c("subt = FDR < 0.05,|coeff| >= 1.25, >20%+ in either group"), 
		ordered=TRUE, p=0.05, FC = 1.25, perc=0.2, extra.cols = c("Treatment"))
```

```
## [1] "0 of 10 samples were removed due to missing data"
```

```
## Default value being used.
```

```
## [1] "Dog will be modeled as a binary categorical predictor variable"
```

```
## Default value being used.
```

```
## it= 0, nll=18.27, log10(eps+1)=Inf, stillActive=135
## it= 1, nll=19.39, log10(eps+1)=Inf, stillActive=8
## it= 2, nll=19.42, log10(eps+1)=Inf, stillActive=3
## it= 3, nll=19.46, log10(eps+1)=Inf, stillActive=2
## it= 4, nll=19.48, log10(eps+1)=Inf, stillActive=1
## it= 5, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it= 6, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it= 7, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it= 8, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it= 9, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=10, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=11, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=12, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=13, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=14, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=15, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=16, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=17, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=18, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=19, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=20, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=21, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=22, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=23, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=24, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=25, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=26, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=27, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=28, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=29, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=30, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=31, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=32, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=33, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=34, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=35, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=36, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=37, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=38, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=39, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=40, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=41, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=42, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=43, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=44, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=45, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=46, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=47, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=48, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=49, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=50, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=51, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=52, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=53, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=54, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=55, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=56, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=57, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=58, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=59, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=60, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=61, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=62, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=63, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=64, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=65, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=66, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=67, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=68, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=69, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=70, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=71, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=72, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=73, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=74, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=75, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=76, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=77, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=78, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=79, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=80, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=81, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=82, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=83, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=84, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=85, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=86, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=87, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=88, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=89, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=90, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=91, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=92, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=93, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=94, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=95, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=96, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=97, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=98, nll=19.49, log10(eps+1)=Inf, stillActive=1
## it=99, nll=19.49, log10(eps+1)=Inf, stillActive=1
## There were  50 OTUs significantly different between B vs. G that met 
##  threshold criteria of p 0.05 absolute FC 1.25 and percentage presence in at least one group of 20 % 
## [1] "writing results and model to file"
## [1] "/scratch/DB/bio/training/16SrRNA/16SrRNA-hex-tutorial/results/downstream_analyses/1_25FC_0.2_Dog_GvsB_OTUs_tax_annot.pdf"
```

![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-22-1.png)![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-22-2.png)

```
## [1] "making heatmap of results"
```

```r
b
```

```
##         percent_positive_group0 percent_positive_group1
## OTU_40                      100                     100
## OTU_43                      100                     100
## OTU_11                      100                     100
## OTU_61                      100                      80
## OTU_87                      100                      20
## OTU_111                      80                       0
## OTU_19                      100                     100
## OTU_77                      100                      80
## OTU_73                      100                     100
## OTU_53                      100                     100
## OTU_18                      100                     100
## OTU_21                      100                     100
## OTU_29                      100                     100
## OTU_33                      100                     100
## OTU_49                      100                     100
## OTU_10                      100                     100
## OTU_46                      100                     100
## OTU_12                      100                     100
## OTU_51                      100                     100
## OTU_59                      100                     100
## OTU_96                       80                     100
## OTU_103                      60                     100
## OTU_150                     100                     100
## OTU_101                     100                     100
## OTU_24                      100                     100
## OTU_140                     100                     100
## OTU_97                      100                     100
## OTU_76                       80                     100
## OTU_54                       80                     100
## OTU_179                      80                     100
## OTU_100                       0                     100
## OTU_144                      40                      20
## OTU_121                     100                     100
## OTU_27                      100                     100
## OTU_8                       100                     100
## OTU_74                       60                     100
## OTU_82                       60                     100
## OTU_99                       40                      40
## OTU_26                      100                     100
## OTU_58                       20                     100
## OTU_69                       60                     100
## OTU_52                       40                     100
## OTU_68                       40                     100
## OTU_72                        0                     100
## OTU_57                       60                     100
## OTU_81                       20                     100
## OTU_36                       40                     100
## OTU_55                       20                     100
## OTU_107                     100                     100
## OTU_42                       40                     100
##         +samples in group 0 +samples in group 1 mean_positive_group0
## OTU_40                    5                   5                  531
## OTU_43                    5                   5                  520
## OTU_11                    5                   5                 3645
## OTU_61                    5                   4                  136
## OTU_87                    5                   1                   36
## OTU_111                   4                   0                   27
## OTU_19                    5                   5                 1142
## OTU_77                    5                   4                   66
## OTU_73                    5                   5                   74
## OTU_53                    5                   5                  152
## OTU_18                    5                   5                 1320
## OTU_21                    5                   5                  582
## OTU_29                    5                   5                 1147
## OTU_33                    5                   5                  103
## OTU_49                    5                   5                  130
## OTU_10                    5                   5                 3376
## OTU_46                    5                   5                 1520
## OTU_12                    5                   5                 1821
## OTU_51                    5                   5                  583
## OTU_59                    5                   5                  105
## OTU_96                    4                   5                    4
## OTU_103                   3                   5                    3
## OTU_150                   5                   5                   19
## OTU_101                   5                   5                   19
## OTU_24                    5                   5                   99
## OTU_140                   5                   5                   32
## OTU_97                    5                   5                    6
## OTU_76                    4                   5                    6
## OTU_54                    4                   5                   32
## OTU_179                   4                   5                    5
## OTU_100                   0                   5                  NaN
## OTU_144                   2                   1                    2
## OTU_121                   5                   5                   14
## OTU_27                    5                   5                  234
## OTU_8                     5                   5                 1510
## OTU_74                    3                   5                    7
## OTU_82                    3                   5                    3
## OTU_99                    2                   2                   22
## OTU_26                    5                   5                   20
## OTU_58                    1                   5                    3
## OTU_69                    3                   5                    4
## OTU_52                    2                   5                    2
## OTU_68                    2                   5                    2
## OTU_72                    0                   5                  NaN
## OTU_57                    3                   5                    2
## OTU_81                    1                   5                    1
## OTU_36                    2                   5                   10
## OTU_55                    1                   5                    4
## OTU_107                   5                   5                    4
## OTU_42                    2                   5                    4
##         mean_positive_group1 oddsRatio      lower       upper     fisherP
## OTU_40                     9  0.000000 0.00000000         Inf 1.000000000
## OTU_43                    13  0.000000 0.00000000         Inf 1.000000000
## OTU_11                    70  0.000000 0.00000000         Inf 1.000000000
## OTU_61                     4       Inf 0.02564066         Inf 1.000000000
## OTU_87                     2       Inf 1.02482226         Inf 0.047619048
## OTU_111                  NaN       Inf 1.02482226         Inf 0.047619048
## OTU_19                    85  0.000000 0.00000000         Inf 1.000000000
## OTU_77                     6       Inf 0.02564066         Inf 1.000000000
## OTU_73                     8  0.000000 0.00000000         Inf 1.000000000
## OTU_53                    17  0.000000 0.00000000         Inf 1.000000000
## OTU_18                   115  0.000000 0.00000000         Inf 1.000000000
## OTU_21                   245  0.000000 0.00000000         Inf 1.000000000
## OTU_29                   241  0.000000 0.00000000         Inf 1.000000000
## OTU_33                    77  0.000000 0.00000000         Inf 1.000000000
## OTU_49                    50  0.000000 0.00000000         Inf 1.000000000
## OTU_10                   804  0.000000 0.00000000         Inf 1.000000000
## OTU_46                   554  0.000000 0.00000000         Inf 1.000000000
## OTU_12                   588  0.000000 0.00000000         Inf 1.000000000
## OTU_51                   149  0.000000 0.00000000         Inf 1.000000000
## OTU_59                    32  0.000000 0.00000000         Inf 1.000000000
## OTU_96                    26  0.000000 0.00000000  39.0005500 1.000000000
## OTU_103                   16  0.000000 0.00000000   5.1183766 0.444444444
## OTU_150                   70  0.000000 0.00000000         Inf 1.000000000
## OTU_101                   86  0.000000 0.00000000         Inf 1.000000000
## OTU_24                  1032  0.000000 0.00000000         Inf 1.000000000
## OTU_140                  147  0.000000 0.00000000         Inf 1.000000000
## OTU_97                    32  0.000000 0.00000000         Inf 1.000000000
## OTU_76                    48  0.000000 0.00000000  39.0005500 1.000000000
## OTU_54                   194  0.000000 0.00000000  39.0005500 1.000000000
## OTU_179                   20  0.000000 0.00000000  39.0005500 1.000000000
## OTU_100                    8  0.000000 0.00000000   0.4353226 0.007936508
## OTU_144                   14  2.414224 0.08474680 195.6529809 1.000000000
## OTU_121                  151  0.000000 0.00000000         Inf 1.000000000
## OTU_27                  1452  0.000000 0.00000000         Inf 1.000000000
## OTU_8                  15389  0.000000 0.00000000         Inf 1.000000000
## OTU_74                    73  0.000000 0.00000000   5.1183766 0.444444444
## OTU_82                    45  0.000000 0.00000000   5.1183766 0.444444444
## OTU_99                    19  1.000000 0.04224561  23.6710987 1.000000000
## OTU_26                   456  0.000000 0.00000000         Inf 1.000000000
## OTU_58                    93  0.000000 0.00000000   0.9757790 0.047619048
## OTU_69                    88  0.000000 0.00000000   5.1183766 0.444444444
## OTU_52                    88  0.000000 0.00000000   2.0268713 0.166666667
## OTU_68                   138  0.000000 0.00000000   2.0268713 0.166666667
## OTU_72                    52  0.000000 0.00000000   0.4353226 0.007936508
## OTU_57                   104  0.000000 0.00000000   5.1183766 0.444444444
## OTU_81                    52  0.000000 0.00000000   0.9757790 0.047619048
## OTU_36                   782  0.000000 0.00000000   2.0268713 0.166666667
## OTU_55                   264  0.000000 0.00000000   0.9757790 0.047619048
## OTU_107                  390  0.000000 0.00000000         Inf 1.000000000
## OTU_42                   617  0.000000 0.00000000   2.0268713 0.166666667
##         fisherAdjP     coeff      pvalues   adjPvalues  Kingdom
## OTU_40   1.0000000 -6.373510 1.161859e-04 1.946113e-03 Bacteria
## OTU_43   1.0000000 -6.116606 3.553180e-08 4.761262e-06 Bacteria
## OTU_11   1.0000000 -4.933445 4.446080e-03 1.715166e-02 Bacteria
## OTU_61   1.0000000 -4.218071 3.696573e-04 2.913770e-03 Bacteria
## OTU_87   0.4945055 -3.748359 5.548120e-04 3.912885e-03 Bacteria
## OTU_111  0.4945055 -3.174336 8.409156e-03 2.691101e-02 Bacteria
## OTU_19   1.0000000 -3.172581 6.922310e-04 4.417093e-03 Bacteria
## OTU_77   1.0000000 -3.142696 1.098060e-02 3.336559e-02 Bacteria
## OTU_73   1.0000000 -3.135823 2.036489e-03 1.091558e-02 Bacteria
## OTU_53   1.0000000 -3.119165 1.374650e-03 8.008833e-03 Bacteria
## OTU_18   1.0000000 -3.088471 4.194520e-05 1.107629e-03 Bacteria
## OTU_21   1.0000000 -2.532059 5.009580e-03 1.864677e-02 Bacteria
## OTU_29   1.0000000 -2.525171 1.153153e-03 7.023750e-03 Bacteria
## OTU_33   1.0000000 -2.470870 8.434795e-03 2.691101e-02 Bacteria
## OTU_49   1.0000000 -1.971867 6.883067e-04 4.417093e-03 Bacteria
## OTU_10   1.0000000 -1.969584 7.249005e-03 2.428417e-02 Bacteria
## OTU_46   1.0000000 -1.849698 1.797311e-03 1.003499e-02 Bacteria
## OTU_12   1.0000000 -1.765926 4.479910e-03 1.715166e-02 Bacteria
## OTU_51   1.0000000 -1.737633 1.286803e-02 3.748513e-02 Bacteria
## OTU_59   1.0000000 -1.694662 3.579636e-03 1.547327e-02 Bacteria
## OTU_96   1.0000000  1.782623 3.303245e-03 1.526327e-02 Bacteria
## OTU_103  1.0000000  1.875474 1.630828e-02 4.379859e-02 Bacteria
## OTU_150  1.0000000  1.928554 1.360546e-02 3.798190e-02 Bacteria
## OTU_101  1.0000000  2.051487 9.561004e-03 2.979476e-02 Bacteria
## OTU_24   1.0000000  2.141545 5.876923e-03 2.118041e-02 Bacteria
## OTU_140  1.0000000  2.222339 6.006384e-03 2.118041e-02 Bacteria
## OTU_97   1.0000000  2.235028 3.419388e-03 1.527327e-02 Bacteria
## OTU_76   1.0000000  2.551360 2.190170e-03 1.128780e-02 Bacteria
## OTU_54   1.0000000  2.784766 7.190695e-03 2.428417e-02 Bacteria
## OTU_179  1.0000000  2.856890 2.538883e-03 1.246244e-02 Bacteria
## OTU_100  0.2678571  2.883321 4.244577e-03 1.715166e-02 Bacteria
## OTU_144  1.0000000  3.083691 1.634276e-02 4.379859e-02 Bacteria
## OTU_121  1.0000000  3.152170 5.164026e-04 3.844331e-03 Bacteria
## OTU_27   1.0000000  3.296192 1.351744e-02 3.798190e-02 Bacteria
## OTU_8    1.0000000  3.446926 2.604093e-03 1.246244e-02 Bacteria
## OTU_74   1.0000000  3.520824 2.496783e-04 2.784229e-03 Bacteria
## OTU_82   1.0000000  3.572533 3.223800e-04 2.913770e-03 Bacteria
## OTU_99   1.0000000  3.927970 1.120486e-02 3.336559e-02 Bacteria
## OTU_26   1.0000000  4.105738 3.395129e-04 2.913770e-03 Bacteria
## OTU_58   0.4945055  4.674182 2.405784e-04 2.784229e-03 Bacteria
## OTU_69   1.0000000  4.879486 3.974350e-03 1.664259e-02 Bacteria
## OTU_52   1.0000000  4.955954 2.132559e-04 2.784229e-03 Bacteria
## OTU_68   1.0000000  4.970360 2.419041e-04 2.784229e-03 Bacteria
## OTU_72   0.2678571  5.046503 2.701118e-04 2.784229e-03 Bacteria
## OTU_57   1.0000000  5.220055 4.959533e-05 1.107629e-03 Bacteria
## OTU_81   0.4945055  5.423223 3.629421e-05 1.107629e-03 Bacteria
## OTU_36   1.0000000  5.803801 6.882302e-05 1.317469e-03 Bacteria
## OTU_55   0.4945055  5.880385 3.425031e-05 1.107629e-03 Bacteria
## OTU_107  1.0000000  6.222838 9.419693e-07 6.311194e-05 Bacteria
## OTU_42   1.0000000  6.234114 3.490600e-04 2.913770e-03 Bacteria
##                 Phylum               Class              Order
## OTU_40      Firmicutes          Clostridia      Clostridiales
## OTU_43      Firmicutes          Clostridia      Clostridiales
## OTU_11  Proteobacteria Gammaproteobacteria  Enterobacteriales
## OTU_61      Firmicutes          Clostridia      Clostridiales
## OTU_87      Firmicutes     Erysipelotrichi Erysipelotrichales
## OTU_111     Firmicutes          Clostridia      Clostridiales
## OTU_19      Firmicutes          Clostridia      Clostridiales
## OTU_77      Firmicutes     Erysipelotrichi Erysipelotrichales
## OTU_73      Firmicutes     Erysipelotrichi Erysipelotrichales
## OTU_53      Firmicutes          Clostridia      Clostridiales
## OTU_18      Firmicutes          Clostridia      Clostridiales
## OTU_21   Bacteroidetes         Bacteroidia      Bacteroidales
## OTU_29      Firmicutes             Bacilli   Turicibacterales
## OTU_33   Bacteroidetes         Bacteroidia      Bacteroidales
## OTU_49      Firmicutes          Clostridia      Clostridiales
## OTU_10      Firmicutes          Clostridia      Clostridiales
## OTU_46      Firmicutes          Clostridia      Clostridiales
## OTU_12      Firmicutes          Clostridia      Clostridiales
## OTU_51      Firmicutes          Clostridia      Clostridiales
## OTU_59      Firmicutes          Clostridia      Clostridiales
## OTU_96      Firmicutes          Clostridia      Clostridiales
## OTU_103     Firmicutes          Clostridia      Clostridiales
## OTU_150     Firmicutes          Clostridia      Clostridiales
## OTU_101   Fusobacteria       Fusobacteriia    Fusobacteriales
## OTU_24      Firmicutes          Clostridia      Clostridiales
## OTU_140   Fusobacteria       Fusobacteriia    Fusobacteriales
## OTU_97      Firmicutes          Clostridia      Clostridiales
## OTU_76      Firmicutes          Clostridia      Clostridiales
## OTU_54   Bacteroidetes         Bacteroidia      Bacteroidales
## OTU_179     Firmicutes          Clostridia      Clostridiales
## OTU_100     Firmicutes     Erysipelotrichi Erysipelotrichales
## OTU_144     Firmicutes          Clostridia      Clostridiales
## OTU_121  Bacteroidetes         Bacteroidia      Bacteroidales
## OTU_27      Firmicutes          Clostridia      Clostridiales
## OTU_8    Bacteroidetes         Bacteroidia      Bacteroidales
## OTU_74      Firmicutes          Clostridia      Clostridiales
## OTU_82      Firmicutes          Clostridia      Clostridiales
## OTU_99      Firmicutes          Clostridia      Clostridiales
## OTU_26   Bacteroidetes         Bacteroidia      Bacteroidales
## OTU_58  Proteobacteria Gammaproteobacteria      Aeromonadales
## OTU_69   Bacteroidetes         Bacteroidia      Bacteroidales
## OTU_52   Bacteroidetes         Bacteroidia      Bacteroidales
## OTU_68  Proteobacteria Gammaproteobacteria      Aeromonadales
## OTU_72  Actinobacteria      Coriobacteriia   Coriobacteriales
## OTU_57  Proteobacteria Gammaproteobacteria      Aeromonadales
## OTU_81      Firmicutes          Clostridia      Clostridiales
## OTU_36   Bacteroidetes         Bacteroidia      Bacteroidales
## OTU_55   Bacteroidetes         Bacteroidia      Bacteroidales
## OTU_107  Bacteroidetes         Bacteroidia      Bacteroidales
## OTU_42   Bacteroidetes         Bacteroidia      Bacteroidales
##                        Family              Genus     Species
## OTU_40        Lachnospiraceae              Dorea            
## OTU_43        Lachnospiraceae              Dorea            
## OTU_11     Enterobacteriaceae                               
## OTU_61  Peptostreptococcaceae                               
## OTU_87    Erysipelotrichaceae                               
## OTU_111       Lachnospiraceae       Epulopiscium            
## OTU_19        Lachnospiraceae     [Ruminococcus]      gnavus
## OTU_77    Erysipelotrichaceae      [Eubacterium]    dolichum
## OTU_73    Erysipelotrichaceae                               
## OTU_53         Clostridiaceae                               
## OTU_18        Lachnospiraceae                               
## OTU_21         Bacteroidaceae        Bacteroides            
## OTU_29      Turicibacteraceae       Turicibacter            
## OTU_33         Bacteroidaceae        Bacteroides            
## OTU_49        Lachnospiraceae                               
## OTU_10        Lachnospiraceae            Blautia    producta
## OTU_46        Lachnospiraceae               <NA>        <NA>
## OTU_12        Ruminococcaceae   Faecalibacterium prausnitzii
## OTU_51                                                      
## OTU_59        Lachnospiraceae     [Ruminococcus]            
## OTU_96        Ruminococcaceae       Oscillospira            
## OTU_103       Ruminococcaceae       Oscillospira            
## OTU_150       Lachnospiraceae                               
## OTU_101      Fusobacteriaceae      Fusobacterium            
## OTU_24        Lachnospiraceae                               
## OTU_140      Fusobacteriaceae      Fusobacterium            
## OTU_97        Ruminococcaceae       Oscillospira            
## OTU_76        Ruminococcaceae       Oscillospira            
## OTU_54         Bacteroidaceae        Bacteroides            
## OTU_179       Lachnospiraceae              Dorea            
## OTU_100   Erysipelotrichaceae                               
## OTU_144       Ruminococcaceae                               
## OTU_121        Prevotellaceae         Prevotella       copri
## OTU_27        Ruminococcaceae       Ruminococcus            
## OTU_8          Prevotellaceae         Prevotella       copri
## OTU_74     [Mogibacteriaceae]                               
## OTU_82        Lachnospiraceae                               
## OTU_99        Veillonellaceae        Megasphaera            
## OTU_26                  S24-7                               
## OTU_58    Succinivibrionaceae Anaerobiospirillum            
## OTU_69         Bacteroidaceae        Bacteroides   uniformis
## OTU_52     Porphyromonadaceae    Parabacteroides            
## OTU_68    Succinivibrionaceae Anaerobiospirillum            
## OTU_72      Coriobacteriaceae      Adlercreutzia            
## OTU_57    Succinivibrionaceae                               
## OTU_81                                                      
## OTU_36   [Paraprevotellaceae]               <NA>        <NA>
## OTU_55                                                      
## OTU_107  [Paraprevotellaceae]       [Prevotella]            
## OTU_42   [Paraprevotellaceae]       [Prevotella]
```

```r
sessionInfo()
```

```
## R version 3.3.0 (2016-05-03)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## Running under: OS X 10.10.5 (Yosemite)
## 
## locale:
## [1] C
## 
## attached base packages:
## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] ROCR_1.0-7           gplots_3.0.1         metagenomeSeq_1.14.2
##  [4] RColorBrewer_1.1-2   glmnet_2.0-5         foreach_1.4.3       
##  [7] Matrix_1.2-6         limma_3.28.5         fifer_1.0           
## [10] MASS_7.3-45          xtable_1.8-2         matrixStats_0.50.2  
## [13] psych_1.6.4          corrplot_0.77        NMF_0.23.4          
## [16] Biobase_2.32.0       BiocGenerics_0.18.0  cluster_2.0.4       
## [19] rngtools_1.2.4       pkgmaker_0.26.6      registry_0.3        
## [22] dplyr_0.5.0          randomForest_4.6-12  vegan_2.4-0         
## [25] lattice_0.20-33      permute_0.9-0        dunn.test_1.3.2     
## [28] gridExtra_2.2.1      ggplot2_2.2.1        phyloseq_1.20.0     
## [31] knitr_1.15.1        
## 
## loaded via a namespace (and not attached):
##  [1] viridis_0.4.0      jsonlite_1.2       viridisLite_0.2.0 
##  [4] splines_3.3.0      gtools_3.5.0       assertthat_0.1    
##  [7] highr_0.6          stats4_3.3.0       robustbase_0.92-6 
## [10] chron_2.3-47       digest_0.6.12      XVector_0.12.0    
## [13] colorspace_1.3-2   plyr_1.8.4         zlibbioc_1.18.0   
## [16] mvtnorm_1.0-5      scales_0.4.1       gdata_2.17.0      
## [19] whisker_0.3-2      tibble_1.3.0       mgcv_1.8-12       
## [22] IRanges_2.6.0      withr_1.0.2        nnet_7.3-12       
## [25] lazyeval_0.2.0     mnormt_1.5-4       survival_2.39-4   
## [28] magrittr_1.5       mclust_5.2         evaluate_0.10     
## [31] doParallel_1.0.10  nlme_3.1-128       class_7.3-14      
## [34] tools_3.3.0        data.table_1.9.6   gridBase_0.4-7    
## [37] trimcluster_0.1-2  stringr_1.2.0      S4Vectors_0.10.1  
## [40] kernlab_0.9-24     munsell_0.4.3      fpc_2.1-10        
## [43] Biostrings_2.40.2  compiler_3.3.0     ade4_1.7-4        
## [46] caTools_1.17.1     rhdf5_2.16.0       grid_3.3.0        
## [49] iterators_1.0.8    biomformat_1.0.2   igraph_1.0.1      
## [52] labeling_0.3       bitops_1.0-6       gtable_0.2.0      
## [55] codetools_0.2-14   multtest_2.28.0    flexmix_2.3-13    
## [58] DBI_0.4-1          reshape2_1.4.2     R6_2.1.2          
## [61] prabclus_2.2-6     KernSmooth_2.23-15 dendextend_1.5.2  
## [64] ape_3.5            modeltools_0.2-21  stringi_1.1.5     
## [67] Rcpp_0.12.10       DEoptimR_1.0-6     diptest_0.75-7
```
