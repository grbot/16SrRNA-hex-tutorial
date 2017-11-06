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
qsub -I -q UCTlong -l walltime=08:00:00
```
Once you are on a compute node you will see that the prompt changes from ```@srvslshpc001``` to ```@srvslshpc60X``` e.g.

```bash
gerrit@srvslshpc001:~> qsub -I -q UCTlong -l walltime=08:00:00
qsub: waiting for job 1598565.srvslshpc001 to start
qsub: job 1598565.srvslshpc001 ready

gerrit@srvslshpc601:~> hostname
srvslshpc601

````

Make an output directory for R downstream analyses
---------------------------------------------------
**NB - replace hpc30 with your user account**

```r
mkdir /researchdata/fhgfs/hpc30/R_downstream
```

```
## Error in eval(expr, envir, enclos): object 'mkdir' not found
```

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
```

```
## Loading required package: ade4
## Loading required package: picante
## Loading required package: ape
## Loading required package: vegan
## Loading required package: permute
## Loading required package: lattice
## This is vegan 2.4-2
## 
## Attaching package: 'vegan'
## 
## The following object is masked from 'package:ade4':
## 
##     cca
## 
## Loading required package: nlme
```

```r
library(ggplot2)
library(gridExtra)
library(dunn.test)
```

```
## Error in library(dunn.test): there is no package called 'dunn.test'
```

```r
library(vegan)
library(randomForest)
```

```
## randomForest 4.6-12
## Type rfNews() to see new features/changes/bug fixes.
## 
## Attaching package: 'randomForest'
## 
## The following object is masked from 'package:gridExtra':
## 
##     combine
```

```r
library(dplyr)
```

```
## Error in library(dplyr): there is no package called 'dplyr'
```
**Import custom functions used in script**


```r
source("/scratch/DB/bio/training/16SrRNA/16SrRNA-hex-tutorial/src/microbiome_custom_functions.R")
```

```
## Warning in grepl("\n", lines, fixed = TRUE): input string 545 is invalid in
## this locale
```

```
## Loading required package: pkgmaker
## Loading required package: registry
## Loading required package: rngtools
## Loading required package: cluster
## NMF - BioConductor layer [OK] | Shared memory capabilities [NO: bigmemory] | Cores 63/64
##   To enable shared memory capabilities, try: install.extras('
## NMF
## ')
## 
## Attaching package: 'NMF'
## 
## The following object is masked from 'package:nlme':
## 
##     coef<-
## 
## The following object is masked from 'package:ape':
## 
##     consensus
```

```
## Error in library(matrixStats): there is no package called 'matrixStats'
```
**Set the working directory and import data**
**NB replace hpc30 with the name that has been given to you in the class**


```r
setwd("/researchdata/fhgfs/hpc30/")
```

```
## Error in setwd("/researchdata/fhgfs/hpc30/"): cannot change working directory
```

```r
inDir <- getwd()#specify input directory
outDir <- paste0(getwd(),"/R_downstream") #specify output directory
phy <- import_biom(BIOMfilename = paste0(inDir,"/otus_table.tax.biom"), 
		verbose = TRUE)#
```

```
## Error in file(con, "r"): cannot open the connection
```

```r
ntaxa(phy) #(number of OTUs)
```

```
## Error in ntaxa(phy): error in evaluating the argument 'physeq' in selecting a method for function 'ntaxa': Error: object 'phy' not found
```

```r
sample_names(phy) <- sub("\\/1","",sample_names(phy))#remove "/1" from filenames
```

```
## Error in sample_names(phy): error in evaluating the argument 'physeq' in selecting a method for function 'sample_names': Error: object 'phy' not found
```

```r
#add phylogenetic tree (.tre file generated in QIIME)
tree <- read_tree_greengenes(paste0(inDir,"/otus_repsetOUT_aligned_pfiltered.tre"))
```

```
## Warning in file(con, "r"): cannot open file '/home/kviljoen/16SrRNA-hex-
## tutorial/downstream/otus_repsetOUT_aligned_pfiltered.tre': No such file or
## directory
```

```
## Error in file(con, "r"): cannot open the connection
```

```r
#merge phy and tree
phy <- merge_phyloseq(phy,tree)
```

```
## Error in merge_phyloseq(phy, tree): object 'phy' not found
```
**Data cleanup**


```r
colnames(tax_table(phy))
```

```
## Error in colnames(tax_table(phy)): error in evaluating the argument 'x' in selecting a method for function 'colnames': Error in tax_table(phy) : 
##   error in evaluating the argument 'object' in selecting a method for function 'tax_table': Error: object 'phy' not found
```

```r
colnames(tax_table(phy)) <-  c("Kingdom", "Phylum" , "Class" , "Order" , "Family" , "Genus", "Species")#e.g. replace "Rank1" with "Kingdom"
```

```
## Error in colnames(tax_table(phy)) <- c("Kingdom", "Phylum", "Class", "Order", : object 'phy' not found
```

```r
#clean taxonomic annotations, at the moment they are for example 'k__Bacteria'; 'p_Firmicutes' - remove k__ and p__ ...
tax_table(phy)[,"Kingdom"] <- sub("k__","",tax_table(phy)[,"Kingdom"])
```

```
## Error in tax_table(phy): error in evaluating the argument 'object' in selecting a method for function 'tax_table': Error: object 'phy' not found
```

```r
tax_table(phy)[,"Phylum"] <- sub("p__","",tax_table(phy)[,"Phylum"])
```

```
## Error in tax_table(phy): error in evaluating the argument 'object' in selecting a method for function 'tax_table': Error: object 'phy' not found
```

```r
tax_table(phy)[,"Class"] <- sub("c__","",tax_table(phy)[,"Class"])
```

```
## Error in tax_table(phy): error in evaluating the argument 'object' in selecting a method for function 'tax_table': Error: object 'phy' not found
```

```r
tax_table(phy)[,"Order"] <- sub("o__","",tax_table(phy)[,"Order"])
```

```
## Error in tax_table(phy): error in evaluating the argument 'object' in selecting a method for function 'tax_table': Error: object 'phy' not found
```

```r
tax_table(phy)[,"Family"] <- sub("f__","",tax_table(phy)[,"Family"])
```

```
## Error in tax_table(phy): error in evaluating the argument 'object' in selecting a method for function 'tax_table': Error: object 'phy' not found
```

```r
tax_table(phy)[,"Genus"] <- sub("g__","",tax_table(phy)[,"Genus"])
```

```
## Error in tax_table(phy): error in evaluating the argument 'object' in selecting a method for function 'tax_table': Error: object 'phy' not found
```

```r
tax_table(phy)[,"Species"] <- sub("s__","",tax_table(phy)[,"Species"])
```

```
## Error in tax_table(phy): error in evaluating the argument 'object' in selecting a method for function 'tax_table': Error: object 'phy' not found
```

```r
t= which(is.na(tax_table(phy)[,"Phylum"])) 
```

```
## Error in tax_table(phy): error in evaluating the argument 'object' in selecting a method for function 'tax_table': Error: object 'phy' not found
```

```r
tax_table(phy) = tax_table(phy)[-t,] #remove rows that don't at least have Phylum-level annotation
```

```
## Error in tax_table(phy): error in evaluating the argument 'object' in selecting a method for function 'tax_table': Error: object 'phy' not found
```
**Import metadata and merge with phyloseq object**


```r
meta <-  read.table(paste0(inDir,"/practice.dataset1.metadata.tsv"), sep = "\t", header =TRUE, row.names=1)
```

```
## Warning in file(file, "rt"): cannot open file '/home/kviljoen/16SrRNA-
## hex-tutorial/downstream/practice.dataset1.metadata.tsv': No such file or
## directory
```

```
## Error in file(file, "rt"): cannot open the connection
```

```r
head(meta)
```

```
## Error in head(meta): object 'meta' not found
```

```r
rownames(meta)
```

```
## Error in rownames(meta): error in evaluating the argument 'x' in selecting a method for function 'rownames': Error: object 'meta' not found
```

```r
head(sample_names(phy))
```

```
## Error in sample_names(phy): error in evaluating the argument 'physeq' in selecting a method for function 'sample_names': Error: object 'phy' not found
```

```r
length(sample_names(phy))#15
```

```
## Error in sample_names(phy): error in evaluating the argument 'physeq' in selecting a method for function 'sample_names': Error: object 'phy' not found
```

```r
length(rownames(meta))#15 (check if same number of samples in .biom file and metadatafile)
```

```
## Error in rownames(meta): error in evaluating the argument 'x' in selecting a method for function 'rownames': Error: object 'meta' not found
```

```r
length(intersect(rownames(meta),sample_names(phy)))#15 (check that the sample names match in all cases)
```

```
## Error in intersect(rownames(meta), sample_names(phy)): error in evaluating the argument 'x' in selecting a method for function 'intersect': Error in rownames(meta) : 
##   error in evaluating the argument 'x' in selecting a method for function 'rownames': Error: object 'meta' not found
```

```r
sample_data(phy) <- meta#assign the metadata to the phyloseq object 'phy' (phyloseq will put these in the right order)
```

```
## Error in eval(expr, envir, enclos): object 'meta' not found
```

```r
nsamples(phy)
```

```
## Error in nsamples(phy): error in evaluating the argument 'physeq' in selecting a method for function 'nsamples': Error: object 'phy' not found
```

```r
str(sample_data(phy))#need to change treatment column to factor variable
```

```
## Error in sample_data(phy): error in evaluating the argument 'object' in selecting a method for function 'sample_data': Error: object 'phy' not found
```

```r
sample_data(phy)[,"Treatment"] <- as.factor(unlist(sample_data(phy)[,"Treatment"]))
```

```
## Error in unlist(sample_data(phy)[, "Treatment"]): error in evaluating the argument 'x' in selecting a method for function 'unlist': Error in sample_data(phy) : 
##   error in evaluating the argument 'object' in selecting a method for function 'sample_data': Error: object 'phy' not found
```
**Save phyloseq object as an .RData file**


```r
save(phy, file = paste0(outDir,"/CBIO_16s_cert.RData")) #Save annotated object as a .RData object for quick reload if required at a later stage
```

```
## Error in save(phy, file = paste0(outDir, "/CBIO_16s_cert.RData")): object 'phy' not found
```

```r
#load(paste0(outDir,"/CBIO_16s_cert.RData")) #this is how you would reload the .RData object 'phy'
```
Explore number of reads per sample, make rarefaction curves and filter data as necessary
-------------------------------------------
**Explore number of reads per sample**

```r
reads <- sample_sums(phy)
```

```
## Error in otu_table(x): error in evaluating the argument 'object' in selecting a method for function 'otu_table': Error: object 'phy' not found
```

```r
length(which(reads<5000))
```

```
## Error in which(reads < 5000): object 'reads' not found
```

```r
raremax <- min(reads)
```

```
## Error in eval(expr, envir, enclos): object 'reads' not found
```

```r
raremax
```

```
## Error in eval(expr, envir, enclos): object 'raremax' not found
```

```r
r=rarecurve(t(otu_table(phy)), step = 100, sample = raremax,xlab = "number of reads/sample", ylab = "number of OTUs",
		label = FALSE, xlim = c(0,100000))
```

```
## Error in t(otu_table(phy)): error in evaluating the argument 'x' in selecting a method for function 't': Error in otu_table(phy) : 
##   error in evaluating the argument 'object' in selecting a method for function 'otu_table': Error: object 'phy' not found
```

```r
pdf(paste0(outDir,"/rarefaction_curve.pdf"))
```

```
## Error in pdf(paste0(outDir, "/rarefaction_curve.pdf")): cannot open file '/home/kviljoen/16SrRNA-hex-tutorial/downstream/R_downstream/rarefaction_curve.pdf'
```

```r
r
```

```
## Error in eval(expr, envir, enclos): object 'r' not found
```

```r
dev.off()
```

```
## null device 
##           1
```

All samples have sufficient sequencing depth for inclusion in downstream analyses. The vertical line in the above plot indicates the sample with the lowest number of reads. Now we will scale data to account for differences in the number of reads/sample and filter rare OTUs that are not of biological interest for the purpose of this analysis (e.g. occurs only in one sample).

**Standardize abundances to median sequence depth**

```r
total = median(sample_sums(phy))
```

```
## Error in otu_table(x): error in evaluating the argument 'object' in selecting a method for function 'otu_table': Error: object 'phy' not found
```

```r
standf = function(x, t=total) round(t * (x / sum(x)))
M.std = transform_sample_counts(phy, standf)
```

```
## Error in fun(1:10): object 'total' not found
```
**Apply mild OTU filter**

Select OTUs where the rowsum for that OTU has at least 20% of samples with a count of 10 each OR where that OTU > 0.001% of the total median count (for cases where the minority of samples may have high counts of a rare OTU)

```r
M.f = filter_taxa(M.std,function(x) sum(x > 10) > (0.02*length(x)) | sum(x) > 0.001*total, TRUE)
```

```
## Error in inherits(x, get.component.classes()): object 'M.std' not found
```

```r
ntaxa(M.f)
```

```
## Error in ntaxa(M.f): error in evaluating the argument 'physeq' in selecting a method for function 'ntaxa': Error: object 'M.f' not found
```
**Basic exploratory plots: alpha- and beta-diversity, barplots, heatmap**
-------------------------------------------
**Alpha diversity by dog**


```r
p <- plot_richness(M.std,x = "Dog",color = "Treatment",measures=c("Shannon"), 
		title = paste0("Standardized to total reads, N=",nsamples(M.std)))+theme(axis.text=element_text(size=16, face="bold"),
				axis.title=element_text(size=16,face="bold"))+geom_point(size=5)
```

```
## Error in otu_table(physeq): error in evaluating the argument 'object' in selecting a method for function 'otu_table': Error: object 'M.std' not found
```

```r
pdf(paste0(outDir,"/alpha_diversity_by_dog_treatment.pdf"))
```

```
## Error in pdf(paste0(outDir, "/alpha_diversity_by_dog_treatment.pdf")): cannot open file '/home/kviljoen/16SrRNA-hex-tutorial/downstream/R_downstream/alpha_diversity_by_dog_treatment.pdf'
```

```r
p
```

```
## Error in eval(expr, envir, enclos): object 'p' not found
```

```r
dev.off()
```

```
## null device 
##           1
```
Is there a significant difference in alpha diversity between dogs irrespective of treatment?

```r
est <- estimate_richness(M.f, split = TRUE, measures = c("Shannon"))
```

```
## Error in otu_table(physeq): error in evaluating the argument 'object' in selecting a method for function 'otu_table': Error: object 'M.f' not found
```

```r
temp <- cbind(est,sample_data(M.f)[,"Dog"])
```

```
## Error in eval(expr, envir, enclos): object 'est' not found
```

```r
head(temp)
```

```
## Error in head(temp): object 'temp' not found
```

```r
t <- kruskal.test(temp[,1]~temp[,2])
```

```
## Error in eval(expr, envir, enclos): object 'temp' not found
```

```r
t
```

```
## standardGeneric for "t" defined from package "base"
## 
## function (x) 
## standardGeneric("t")
## <environment: 0x5e29b18>
## Methods may be defined for arguments: x
## Use  showMethods("t")  for currently available ones.
```

```r
dunn.test(temp[,1],temp[,2])#post-hoc testing to see which dogs are different
```

```
## Error in eval(expr, envir, enclos): could not find function "dunn.test"
```
Dog G has higher alpha diversity than dogs K and B irrespective of treatment, but this difference is not significant

**Alpha diversity by treatment**

```r
p <- plot_richness(M.std,x = "Treatment",color = "Dog",measures=c("Shannon"), 
				title = paste0("Standardized to total reads, N=",nsamples(M.std)))+theme(axis.text=element_text(size=16, face="bold"),
				axis.title=element_text(size=16,face="bold"))+geom_point(size=5)
```

```
## Error in otu_table(physeq): error in evaluating the argument 'object' in selecting a method for function 'otu_table': Error: object 'M.std' not found
```

```r
pdf(paste0(outDir,"/alpha_diversity_by_treatment_dog.pdf"))
```

```
## Error in pdf(paste0(outDir, "/alpha_diversity_by_treatment_dog.pdf")): cannot open file '/home/kviljoen/16SrRNA-hex-tutorial/downstream/R_downstream/alpha_diversity_by_treatment_dog.pdf'
```

```r
p
```

```
## Error in eval(expr, envir, enclos): object 'p' not found
```

```r
dev.off()
```

```
## null device 
##           1
```
Are there significant differences in alpha diversity by treatment?

```r
temp <- cbind(est,sample_data(M.f)[,"Treatment"])
```

```
## Error in eval(expr, envir, enclos): object 'est' not found
```

```r
head(temp)
```

```
## Error in head(temp): object 'temp' not found
```

```r
t <- kruskal.test(temp[,1]~temp[,2])
```

```
## Error in eval(expr, envir, enclos): object 'temp' not found
```

```r
t
```

```
## standardGeneric for "t" defined from package "base"
## 
## function (x) 
## standardGeneric("t")
## <environment: 0x5e29b18>
## Methods may be defined for arguments: x
## Use  showMethods("t")  for currently available ones.
```

```r
dunn.test(temp[,1],temp[,2])
```

```
## Error in eval(expr, envir, enclos): could not find function "dunn.test"
```

**Beta diversity using NMDS with Bray-Curtis as distance metric**

```r
set.seed(2)
GP.ord.BC <- ordinate(M.f, "NMDS", "bray", k=2, trymax=100)
```

```
## Error in ordinate(M.f, "NMDS", "bray", k = 2, trymax = 100): object 'M.f' not found
```

```r
GP.ord.BC
```

```
## Error in eval(expr, envir, enclos): object 'GP.ord.BC' not found
```

```r
color = c("Treatment")
shape = c("Dog")
title=c("NMDS of 16S microbiome,Bray-Curtis distance,k=2")
MDS = plot_ordination(M.f, GP.ord.BC, color = color,shape=shape, 
		title = title)
```

```
## Error in plot_ordination(M.f, GP.ord.BC, color = color, shape = shape, : object 'M.f' not found
```

```r
MDS.1  = MDS +theme(axis.text=element_text(size=16, face="bold"),
				axis.title=element_text(size=18,face="bold"), legend.title=element_text(size=14))+
		theme_bw()+labs(color=color, shape=shape)+geom_point(size=5)
```

```
## Error in eval(expr, envir, enclos): object 'MDS' not found
```

```r
pdf(paste0(outDir,"/NMDS_Dogs_treatment_Bray_Curtis.pdf"),8,5)
```

```
## Error in pdf(paste0(outDir, "/NMDS_Dogs_treatment_Bray_Curtis.pdf"), 8, : cannot open file '/home/kviljoen/16SrRNA-hex-tutorial/downstream/R_downstream/NMDS_Dogs_treatment_Bray_Curtis.pdf'
```

```r
MDS.1
```

```
## Error in eval(expr, envir, enclos): object 'MDS.1' not found
```

```r
dev.off()
```

```
## null device 
##           1
```
**Beta diversity using NMDS with Unifrac as distance metric**

```r
GP.ord.U <- ordinate(M.f, "NMDS", "unifrac")
```

```
## Error in ordinate(M.f, "NMDS", "unifrac"): object 'M.f' not found
```

```r
GP.ord.U
```

```
## Error in eval(expr, envir, enclos): object 'GP.ord.U' not found
```

```r
color = c("Treatment")
shape = c("Dog")

title=c("NMDS of 16S microbiome, Unifrac distance, k=2")

MDS = plot_ordination(M.f, GP.ord.U, color = color, shape=shape, 
		title = title)
```

```
## Error in plot_ordination(M.f, GP.ord.U, color = color, shape = shape, : object 'M.f' not found
```

```r
MDS.1  = MDS +theme(axis.text=element_text(size=16, face="bold"),
				axis.title=element_text(size=18,face="bold"), legend.title=element_text(size=14))+
		theme_bw()+labs(color=color)+geom_point(size=5)
```

```
## Error in eval(expr, envir, enclos): object 'MDS' not found
```

```r
pdf(paste0(outDir,"/NMDS_Dogs_treatment_Unifrac.pdf"),8,5)
```

```
## Error in pdf(paste0(outDir, "/NMDS_Dogs_treatment_Unifrac.pdf"), 8, 5): cannot open file '/home/kviljoen/16SrRNA-hex-tutorial/downstream/R_downstream/NMDS_Dogs_treatment_Unifrac.pdf'
```

```r
MDS.1
```

```
## Error in eval(expr, envir, enclos): object 'MDS.1' not found
```

```r
dev.off()
```

```
## null device 
##           1
```
**Create a heatmap of taxa merged at the lowest available taxonomic level**

```r
M.phy <- tax_glom.kv(M.f)#this function is available in the 'microbiome_custom_functions.R' script loaded at the beginning of this script
```

```
## Error in eval(expr, envir, enclos): could not find function "tax_glom.kv"
```

```r
ntaxa(M.phy)
```

```
## Error in ntaxa(M.phy): error in evaluating the argument 'physeq' in selecting a method for function 'ntaxa': Error: object 'M.phy' not found
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
```

```
## Error in distance(M.phy, method = "bray", type = "samples"): object 'M.phy' not found
```

```r
clust.res<-hclust(diss)
```

```
## Error in hclust(diss): object 'diss' not found
```

```r
sample.order = clust.res$order
```

```
## Error in eval(expr, envir, enclos): object 'clust.res' not found
```

```r
#heatmap is output to file (the heatmap.k function can be found in the 'microbiome_custom_functions.R' script)
hm = heatmap.k(physeq= M.phy,
		annot.cols = c(1,2),
		main = main,filename = f,colours=colours,Colv = sample.order,labrow = TRUE)	
```

```
## Error in eval(expr, envir, enclos): could not find function "heatmap.k"
```

```r
print(hm)
```

```
## Error in print(hm): object 'hm' not found
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
```

```
## Error in eval(expr, envir, enclos): could not find function "bar.plots"
```

```r
print(barplot)
```

```
## function (height, ...) 
## UseMethod("barplot")
## <bytecode: 0xf933b8>
## <environment: namespace:graphics>
```

Detect taxa/OTUs that differ significantly by Dog
-------------------------------------------
convert phyloseq object to metagenomeSeq obj. NB use raw data not standardized:

```r
Mraw.f = filter_taxa(phy,function(x) sum(x > 10) > (0.02*length(x)) | sum(x) > 0.001*total, TRUE)
```

```
## Error in inherits(x, get.component.classes()): object 'phy' not found
```

```r
ntaxa(Mraw.f)
```

```
## Error in ntaxa(Mraw.f): error in evaluating the argument 'physeq' in selecting a method for function 'ntaxa': Error: object 'Mraw.f' not found
```

```r
MGS=make_metagenomeSeq(Mraw.f)
```

```
## Error in eval(expr, envir, enclos): could not find function "make_metagenomeSeq"
```

```r
MGS
```

```
## Error in eval(expr, envir, enclos): object 'MGS' not found
```

**Use Random forests analysis to detect taxa that are good predictors of Dog**

Example used: Dog G vs. Dog B (all treatment points)


```r
sub.index <- sample_names(M.f)[sample_data(M.f)[,"Dog"] != "K"]
```

```
## Error in sample_names(M.f): error in evaluating the argument 'physeq' in selecting a method for function 'sample_names': Error: object 'M.f' not found
```

```r
phy.temp <- prune_samples(sub.index, M.f)
```

```
## Error in prune_samples(sub.index, M.f): error in evaluating the argument 'samples' in selecting a method for function 'prune_samples': Error: object 'sub.index' not found
```

```r
nsamples(phy.temp)
```

```
## Error in nsamples(phy.temp): error in evaluating the argument 'physeq' in selecting a method for function 'nsamples': Error: object 'phy.temp' not found
```

```r
RF.k(data = phy.temp, var = "Dog", ntree=10000, cv.fold=10, outDir = outDir, Nfeatures.validation = 3)
```

```
## Error in eval(expr, envir, enclos): could not find function "RF.k"
```

The class error rates are 0% (even one OTU enough to discriminate between Dog G and B?)

What if we used merged OTUs?

```r
merged.phy <- tax_glom.kv(phy.temp)
```

```
## Error in eval(expr, envir, enclos): could not find function "tax_glom.kv"
```

```r
RF.k(data = merged.phy, var = "Dog", ntree=10000, cv.fold=10, outDir = outDir, Nfeatures.validation = 3, descriptor = "merged_OTUs") #for details on RF.k() see microbiome_custom_functions.R file
```

```
## Error in eval(expr, envir, enclos): could not find function "RF.k"
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
## Error in eval(expr, envir, enclos): could not find function "super.fitZig.kv"
```

```r
print(a)
```

```
## Error in print(a): object 'a' not found
```
Now again compare dog G vs. dog B with differential abundance testing (individual OTUs as opposed to merged this time)


```r
b = super.fitZig.kv(physeq = phy.temp,factor = "Dog",outDir = outDir,FileName =c("1_25FC_0.2_Dog_GvsB_OTUs"),
		heatmap.descriptor=c("tax_annot"), main=c("Dog G vs. B, OTUs"), subt=c("subt = FDR < 0.05,|coeff| >= 1.25, >20%+ in either group"), 
		ordered=TRUE, p=0.05, FC = 1.25, perc=0.2, extra.cols = c("Treatment"))
```

```
## Error in eval(expr, envir, enclos): could not find function "super.fitZig.kv"
```

```r
b
```

```
## Error in eval(expr, envir, enclos): object 'b' not found
```

```r
sessionInfo()
```

```
## R version 3.0.2 (2013-09-25)
## Platform: x86_64-unknown-linux-gnu (64-bit)
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] NMF_0.20.5          Biobase_2.22.0      BiocGenerics_0.8.0 
##  [4] cluster_2.0.5       rngtools_1.2.4      pkgmaker_0.22      
##  [7] registry_0.2        randomForest_4.6-12 gridExtra_2.2.1    
## [10] ggplot2_1.0.1       phyloseq_1.6.1      picante_1.6-2      
## [13] nlme_3.1-131        vegan_2.4-2         lattice_0.20-34    
## [16] permute_0.9-4       ape_3.0-11          ade4_1.7-5         
## [19] knitr_1.11         
## 
## loaded via a namespace (and not attached):
##  [1] biom_0.3.12        Biostrings_2.30.1  codetools_0.2-15  
##  [4] colorspace_1.3-2   digest_0.6.12      doParallel_1.0.8  
##  [7] evaluate_0.10      foreach_1.4.3      formatR_1.4       
## [10] grid_3.0.2         gridBase_0.4-7     gtable_0.2.0      
## [13] igraph_1.0.1       IRanges_1.20.7     iterators_1.0.8   
## [16] magrittr_1.5       MASS_7.3-29        Matrix_1.2-8      
## [19] mgcv_1.8-17        multtest_2.18.0    munsell_0.4.3     
## [22] plyr_1.8.1         proto_1.0.0        RColorBrewer_1.0-5
## [25] Rcpp_0.12.9        reshape2_1.4.2     RJSONIO_1.3-0     
## [28] scales_0.4.1       splines_3.0.2      stats4_3.0.2      
## [31] stringi_1.1.2      stringr_1.2.0      survival_2.40-1   
## [34] tools_3.0.2        xtable_1.8-2       XVector_0.2.0
```
