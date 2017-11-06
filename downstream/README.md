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
library(ggplot2)
library(gridExtra)
library(dunn.test)
library(vegan)
library(randomForest)
library(dplyr)
```
**Import custom functions used in script**


```r
source("/scratch/DB/bio/training/16SrRNA/16SrRNA-hex-tutorial/src/microbiome_custom_functions.R")
```

```
## Warning in file(filename, "r", encoding = encoding): cannot open
## file '/scratch/DB/bio/training/16SrRNA/16SrRNA-hex-tutorial/src/
## microbiome_custom_functions.R': No such file or directory
```

```
## Error in file(filename, "r", encoding = encoding): cannot open the connection
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
## Error in read_biom(biom_file = BIOMfilename): Both attempts to read input file:
## /Users/katielennard/Documents/Academic/Postdoc/Projects/cbio_16S_pipeline_certification/github_repo/otus_table.tax.biom
## either as JSON (BIOM-v1) or HDF5 (BIOM-v2).
## Check file path, file name, file itself, then try again.
```

```r
ntaxa(phy) #(number of OTUs)
```

```
## Error in ntaxa(phy): object 'phy' not found
```

```r
sample_names(phy) <- sub("\\/1","",sample_names(phy))#remove "/1" from filenames
```

```
## Error in sample_names(phy): object 'phy' not found
```

```r
#add phylogenetic tree (.tre file generated in QIIME)
tree <- read_tree_greengenes(paste0(inDir,"/otus_repsetOUT_aligned_pfiltered.tre"))
```

```
## Warning in file(con, "r"): cannot open file '/Users/katielennard/Documents/
## Academic/Postdoc/Projects/cbio_16S_pipeline_certification/github_repo/
## otus_repsetOUT_aligned_pfiltered.tre': No such file or directory
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
## Error in tax_table(phy): object 'phy' not found
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
## Error in tax_table(phy): object 'phy' not found
```

```r
tax_table(phy)[,"Phylum"] <- sub("p__","",tax_table(phy)[,"Phylum"])
```

```
## Error in tax_table(phy): object 'phy' not found
```

```r
tax_table(phy)[,"Class"] <- sub("c__","",tax_table(phy)[,"Class"])
```

```
## Error in tax_table(phy): object 'phy' not found
```

```r
tax_table(phy)[,"Order"] <- sub("o__","",tax_table(phy)[,"Order"])
```

```
## Error in tax_table(phy): object 'phy' not found
```

```r
tax_table(phy)[,"Family"] <- sub("f__","",tax_table(phy)[,"Family"])
```

```
## Error in tax_table(phy): object 'phy' not found
```

```r
tax_table(phy)[,"Genus"] <- sub("g__","",tax_table(phy)[,"Genus"])
```

```
## Error in tax_table(phy): object 'phy' not found
```

```r
tax_table(phy)[,"Species"] <- sub("s__","",tax_table(phy)[,"Species"])
```

```
## Error in tax_table(phy): object 'phy' not found
```

```r
t= which(is.na(tax_table(phy)[,"Phylum"])) 
```

```
## Error in tax_table(phy): object 'phy' not found
```

```r
tax_table(phy) = tax_table(phy)[-t,] #remove rows that don't at least have Phylum-level annotation
```

```
## Error in tax_table(phy): object 'phy' not found
```
**Import metadata and merge with phyloseq object**


```r
meta <-  read.table(paste0(inDir,"/practice.dataset1.metadata.tsv"), sep = "\t", header =TRUE, row.names=1)
```

```
## Warning in file(file, "rt"): cannot open file '/Users/katielennard/
## Documents/Academic/Postdoc/Projects/cbio_16S_pipeline_certification/
## github_repo/practice.dataset1.metadata.tsv': No such file or directory
```

```
## Error in file(file, "rt"): cannot open the connection
```

```r
head(meta)
```

```
## Sample Data:        [6 samples by 95 sample variables]:
##      CT NG TV MG HSV.1 HSV.2           BV Nugent.score STI.any  pH
## SW16  0  0  0  0     0     0 Intermediate            4       0 6.1
## SW17  0  0  0  0     0     0     Negative            0       0 4.1
## SW18  0  0  0  0     0     0     Positive            7       0  NA
## SW23  0  0  0  0     0     0     Negative            2       0 4.4
## SW24  0  0  0  0     0     0 Intermediate            4       0 4.1
## SW25  0  0  0  0     0     0     Positive           10       0 5.0
##      location     IL.1a   IL.2Ra      IL.3 IL.12p40     IL.16     IL.18
## SW16      JHB 11.393122 9.910942 11.756035 12.74412 14.314116 12.554073
## SW17      JHB  7.326879 6.430453  8.509973 10.84584  7.322379  4.436295
## SW18      JHB 10.140319 6.720415  9.114653 10.81045  7.374170  5.839204
## SW23      JHB 11.992602 9.760470 11.808461 12.06950 11.249854 14.448000
## SW24      JHB        NA       NA        NA       NA        NA        NA
## SW25      JHB  9.701913 9.326205 11.173708 11.89059  9.500941 13.022238
##         CTACK      GROa       HGF   IFN.a2      LIF    MCP.3     M.CSF
## SW16 9.226894 13.057400 15.449374 8.473503 8.968091 7.948367  9.183759
## SW17 6.157852  9.770168  6.066089 5.758889 5.149747 4.169368  8.347178
## SW18 6.455327  8.437752  6.508587 6.104337 5.865424 4.169368  9.396284
## SW23 9.062451 13.057400 13.662073 8.503428 9.051345 7.917074  8.089318
## SW24       NA        NA        NA       NA       NA       NA        NA
## SW25 8.286327 11.479376 11.418011 8.114523 8.340518 7.284477 10.879392
##            MIF      MIG     b.NGF      SCF    SCGF.b   SDF.1a    TNF.b
## SW16 17.575438 14.45616 5.7839804 9.749450 17.312818 9.512839 5.383704
## SW17  6.502235 12.34993 0.2630344 5.586465  5.873205 4.638644 3.842979
## SW18 10.715576 11.90599 0.6322682 5.236493  5.873205 7.410663 4.285402
## SW23 16.323061 11.20460 5.6438562 7.876210 11.807033 9.238763 6.818902
## SW24        NA       NA        NA       NA        NA       NA       NA
## SW25 16.367848 14.18324 4.7114949 7.905989 11.273533 8.839361 6.029011
##          TRAIL     IL.1b   IL.1ra     IL.2      IL.4     IL.5      IL.6
## SW16 13.860292 12.849180 15.55031 6.403438  4.984134 4.596935 10.508686
## SW17  6.034524 -1.321928 13.05880 5.580447 -1.321928 2.880776  2.687520
## SW18  6.430453  8.219411 15.29601 3.622776  2.916477 2.880776  3.446256
## SW23 11.415161 10.641826 15.11836 3.622776  2.981853 2.880776  2.687520
## SW24        NA        NA       NA       NA        NA       NA        NA
## SW25 10.178478  8.145677 15.09247 3.622776  2.446256 2.880776  5.727920
##          IL.7      IL.8     IL.9    IL.10 IL.12.p70.    IL.13    IL.17
## SW16 5.063934 16.937948 7.357112 9.859379   8.995202 6.491853 8.649795
## SW17 4.259272  4.778734 3.482892 6.577429   5.575032 2.307355 8.906740
## SW18 5.063934  9.422170 5.364572 7.506605   9.826310 2.307355 8.906740
## SW23 2.733354 13.655654 6.010108 7.635537   5.575032 2.307355 8.906740
## SW24       NA        NA       NA       NA         NA       NA       NA
## SW25 2.733354 10.130892 5.446256 7.499447   5.575032 2.307355 8.906740
##       Eotaxin FGF.basic     G.CSF   GM.CSF    IFN.g     IP.10     MCP.1
## SW16 7.765203  7.947783 10.527135 6.404290 8.952304 11.610471 11.633018
## SW17 7.364135  7.858913  7.106432 8.557655 6.475733  7.994071  6.292782
## SW18 8.158610  7.591335  8.709945 8.308111 7.365885  9.234578  7.702173
## SW23 7.871289  8.772480  8.543032 8.067703 7.705632 10.082947  6.075747
## SW24       NA        NA        NA       NA       NA        NA        NA
## SW25 7.871289  8.685800 12.309690 8.179412 7.641329 15.209564  9.278333
##         MIP.1a   PDGF.bb    MIP.1b    RANTES     TNF.a      VEGF
## SW16 15.572945 10.153235 12.438610 11.917745 12.989705 12.893083
## SW17  2.153805 -1.234465  4.296457  4.993221  4.274262 12.849806
## SW18  2.608809  8.025416  7.657140  7.967514  5.546894 14.073275
## SW23  4.422906  7.369815  9.121145  6.418696  7.485829  7.350464
## SW24        NA        NA        NA        NA        NA        NA
## SW25  3.892391  7.625344  6.981282  7.001127  6.675957  7.350464
##      contraceptive Oestradiol.nmol.l. Progesterone..nmol.l.
## SW16          DMPA                 NA                    NA
## SW17   Male Condom                 NA                    NA
## SW18  Nur isterate                 NA                    NA
## SW23  Nur isterate             0.0472                   0.5
## SW24  Nur isterate             0.1254                   1.7
## SW25  Nur isterate             0.1436                   1.2
##      Luteinising.Hormone..IU.l. High.Med.Low.PSA CD4.freq CD4.CCR5.freq
## SW16                         NA                L       NA            NA
## SW17                         NA                L       NA            NA
## SW18                         NA                M       NA            NA
## SW23                        6.6                L       NA            NA
## SW24                        6.1                L       NA            NA
## SW25                       11.1                L       NA            NA
##      CD4.CD38.freq CD4.CD38.HLADR.freq CD4.HLADR.freq CD4.Ki67.freq
## SW16            NA                  NA             NA            NA
## SW17            NA                  NA             NA            NA
## SW18            NA                  NA             NA            NA
## SW23            NA                  NA             NA            NA
## SW24            NA                  NA             NA            NA
## SW25            NA                  NA             NA            NA
##      CD8.freq CD8.CCR5.freq CD8.CD38.freq CD8.CD38.HLADR.freq
## SW16       NA            NA            NA                  NA
## SW17       NA            NA            NA                  NA
## SW18       NA            NA            NA                  NA
## SW23       NA            NA            NA                  NA
## SW24       NA            NA            NA                  NA
## SW25       NA            NA            NA                  NA
##      CD8.HLADR.freq CD8.Ki67.freq CD4..CCR5..MFI CD8..CCR5..MFI
## SW16             NA            NA             NA             NA
## SW17             NA            NA             NA             NA
## SW18             NA            NA             NA             NA
## SW23             NA            NA             NA             NA
## SW24             NA            NA             NA             NA
## SW25             NA            NA             NA             NA
##      CD4..CCR5..CD38. CD4..CCR5..CD38.HLADR. CD4..CCR5..HLADR.
## SW16               NA                     NA                NA
## SW17               NA                     NA                NA
## SW18               NA                     NA                NA
## SW23               NA                     NA                NA
## SW24               NA                     NA                NA
## SW25               NA                     NA                NA
##      CD4..CCR5..Ki67. CD8..CCR5..CD38. CD8..CCR5..CD38.HLADR.
## SW16               NA               NA                     NA
## SW17               NA               NA                     NA
## SW18               NA               NA                     NA
## SW23               NA               NA                     NA
## SW24               NA               NA                     NA
## SW25               NA               NA                     NA
##      CD8..CCR5..HLADR. CD8..CCR5..Ki67. HPV.risk    BMI
## SW16                NA               NA      Neg 21.219
## SW17                NA               NA      Neg 22.481
## SW18                NA               NA      Neg     NA
## SW23                NA               NA     High 22.151
## SW24                NA               NA      Low 22.151
## SW25                NA               NA     High 29.372
##                        Anahtar.subtype Ethnicity     injectable Age
## SW16                           L.iners      <NA>     Injectable  18
## SW17 Lactobacillus (excluding L.iners)      <NA> No injectables  19
## SW18                           L.iners      <NA>     Injectable  NA
## SW23                           L.iners      <NA>     Injectable  18
## SW24                           L.iners      <NA>     Injectable  18
## SW25                        Prevotella      <NA>     Injectable  19
##      Functional.subtype Compositional.subtype Inflammation  HC
## SW16                 F2                    C1         High Yes
## SW17                 F1                    C2          Low  No
## SW18                 F1                    C3          Low Yes
## SW23                 F1                    C3         High Yes
## SW24                 F1                    C3         <NA> Yes
## SW25                 F3                    C1         High Yes
```

```r
rownames(meta)
```

```
##   [1] "SW16"  "SW17"  "SW18"  "SW23"  "SW24"  "SW25"  "SW26"  "SW27" 
##   [9] "SW28"  "SW29"  "SW31"  "SW32"  "SW33"  "SW34"  "SW35"  "SW36" 
##  [17] "SW37"  "SW38"  "SW39"  "SW42"  "SW43"  "SW44"  "SW45"  "SW46" 
##  [25] "SW47"  "SW48"  "SW49"  "SW50"  "SW53"  "SW54"  "SW55"  "SW56" 
##  [33] "SW57"  "SW58"  "SW59"  "SW60"  "SW61"  "SW62"  "SW63"  "SW64" 
##  [41] "SW65"  "SW66"  "SW67"  "SW68"  "SW69"  "SW70"  "SW71"  "SW72" 
##  [49] "SW73"  "SW74"  "SW75"  "SW76"  "SW77"  "SW78"  "SW79"  "SW80" 
##  [57] "SW81"  "SW82"  "SW83"  "SW84"  "SW87"  "SW88"  "SW89"  "SW90" 
##  [65] "SW91"  "SW92"  "SW93"  "SW94"  "SW95"  "SW96"  "SW97"  "SW98" 
##  [73] "SW99"  "SW100" "SW101" "SW102" "SW103" "SW104" "W2"    "W4"   
##  [81] "W6"    "W7"    "W8"    "W9"    "W10"   "W11"   "W13"   "W15"  
##  [89] "W17"   "W21"   "W23"   "W24"   "W27"   "W28"   "W30"   "W31"  
##  [97] "W32"   "W33"   "W35"   "W36"   "W37"   "W39"   "W40"   "W43"  
## [105] "W45"   "W48"   "W51"   "W52"   "W53"   "W54"   "W56"   "W57"  
## [113] "W63"   "W65"   "W68"   "W70"   "W71"   "W72"   "W73"   "W74"  
## [121] "W77"   "W79"   "W81"   "W82"   "W84"   "W85"   "W86"   "W87"  
## [129] "W88"   "W89"   "W91"   "W92"   "W94"   "W95"   "W96"   "W97"  
## [137] "W99"   "W100"  "W101"  "W102"  "W106"  "W108"  "W110"  "W112" 
## [145] "W113"  "W114"  "W116"  "W117"  "W119"  "W120"  "W122"  "W124" 
## [153] "W125"  "W126"  "W127"  "W128"  "W129"  "W130"  "W131"  "W132" 
## [161] "W136"  "W137"  "W139"  "W141"  "W147"  "W149"  "W152"  "W154"
```

```r
head(sample_names(phy))
```

```
## Error in sample_names(phy): object 'phy' not found
```

```r
length(sample_names(phy))#15
```

```
## Error in sample_names(phy): object 'phy' not found
```

```r
length(rownames(meta))#15 (check if same number of samples in .biom file and metadatafile)
```

```
## [1] 168
```

```r
length(intersect(rownames(meta),sample_names(phy)))#15 (check that the sample names match in all cases)
```

```
## Error in sample_names(phy): object 'phy' not found
```

```r
sample_data(phy) <- meta#assign the metadata to the phyloseq object 'phy' (phyloseq will put these in the right order)
```

```
## Error in sample_data(phy) <- meta: object 'phy' not found
```

```r
nsamples(phy)
```

```
## Error in nsamples(phy): object 'phy' not found
```

```r
str(sample_data(phy))#need to change treatment column to factor variable
```

```
## Error in sample_data(phy): object 'phy' not found
```

```r
sample_data(phy)[,"Treatment"] <- as.factor(unlist(sample_data(phy)[,"Treatment"]))
```

```
## Error in sample_data(phy): object 'phy' not found
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
## Error in otu_table(x): object 'phy' not found
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
## Error in otu_table(phy): object 'phy' not found
```

```r
pdf(paste0(outDir,"/rarefaction_curve.pdf"))
```

```
## Error in pdf(paste0(outDir, "/rarefaction_curve.pdf")): cannot open file '/Users/katielennard/Documents/Academic/Postdoc/Projects/cbio_16S_pipeline_certification/github_repo/R_downstream/rarefaction_curve.pdf'
```

```r
r
```

```
## function (stringA, stringB, data.names, names = F) 
## {
##     start = which(data.names == stringA)
##     end = which(data.names == stringB)
##     if (length(start) < 1) {
##         warning(paste("\nCould not find the variable", stringA, 
##             "in the names"))
##         stop()
##     }
##     if (length(end) < 1) {
##         warning(paste("\nCould not find the variable", stringB, 
##             "in the names"))
##         stop()
##     }
##     if (names) {
##         return(data.names[start:end])
##     }
##     else {
##         return(start:end)
##     }
## }
## <environment: namespace:fifer>
```

```r
dev.off()
```

```
## RStudioGD 
##         2
```

All samples have sufficient sequencing depth for inclusion in downstream analyses. The vertical line in the above plot indicates the sample with the lowest number of reads. Now we will scale data to account for differences in the number of reads/sample and filter rare OTUs that are not of biological interest for the purpose of this analysis (e.g. occurs only in one sample).

**Standardize abundances to median sequence depth**

```r
total = median(sample_sums(phy))
```

```
## Error in otu_table(x): object 'phy' not found
```

```r
standf = function(x, t=total) round(t * (x / sum(x)))
M.std = transform_sample_counts(phy, standf)
```

```
## Error in taxa_are_rows(physeq): object 'phy' not found
```
**Apply mild OTU filter**

Select OTUs where the rowsum for that OTU has at least 20% of samples with a count of 10 each OR where that OTU > 0.001% of the total median count (for cases where the minority of samples may have high counts of a rare OTU)

```r
M.f = filter_taxa(M.std,function(x) sum(x > 10) > (0.02*length(x)) | sum(x) > 0.001*total, TRUE)
ntaxa(M.f)
```

```
## [1] 199
```
**Basic exploratory plots: alpha- and beta-diversity, barplots, heatmap**
-------------------------------------------
**Alpha diversity by dog**


```r
p <- plot_richness(M.std,x = "Dog",color = "Treatment",measures=c("Shannon"), 
		title = paste0("Standardized to total reads, N=",nsamples(M.std)))+theme(axis.text=element_text(size=16, face="bold"),
				axis.title=element_text(size=16,face="bold"))+geom_point(size=5)

pdf(paste0(outDir,"/alpha_diversity_by_dog_treatment.pdf"))
```

```
## Error in pdf(paste0(outDir, "/alpha_diversity_by_dog_treatment.pdf")): cannot open file '/Users/katielennard/Documents/Academic/Postdoc/Projects/cbio_16S_pipeline_certification/github_repo/R_downstream/alpha_diversity_by_dog_treatment.pdf'
```

```r
p
```

```
## Error in eval(expr, envir, enclos): object 'Treatment' not found
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png)

```r
dev.off()
```

```
## RStudioGD 
##         2
```
Is there a significant difference in alpha diversity between dogs irrespective of treatment?

```r
est <- estimate_richness(M.f, split = TRUE, measures = c("Shannon"))
temp <- cbind(est,sample_data(M.f)[,"Dog"])
```

```
## Error in `[.data.frame`(data.frame(x), i, j, drop = FALSE): undefined columns selected
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
## <environment: 0x10619c4a0>
## Methods may be defined for arguments: x
## Use  showMethods("t")  for currently available ones.
```

```r
dunn.test(temp[,1],temp[,2])#post-hoc testing to see which dogs are different
```

```
## Error in dunn.test(temp[, 1], temp[, 2]): object 'temp' not found
```
Dog G has higher alpha diversity than dogs K and B irrespective of treatment, but this difference is not significant

**Alpha diversity by treatment**

```r
p <- plot_richness(M.std,x = "Treatment",color = "Dog",measures=c("Shannon"), 
				title = paste0("Standardized to total reads, N=",nsamples(M.std)))+theme(axis.text=element_text(size=16, face="bold"),
				axis.title=element_text(size=16,face="bold"))+geom_point(size=5)

pdf(paste0(outDir,"/alpha_diversity_by_treatment_dog.pdf"))
```

```
## Error in pdf(paste0(outDir, "/alpha_diversity_by_treatment_dog.pdf")): cannot open file '/Users/katielennard/Documents/Academic/Postdoc/Projects/cbio_16S_pipeline_certification/github_repo/R_downstream/alpha_diversity_by_treatment_dog.pdf'
```

```r
p
```

```
## Error in eval(expr, envir, enclos): object 'Dog' not found
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13-1.png)

```r
dev.off()
```

```
## RStudioGD 
##         2
```
Are there significant differences in alpha diversity by treatment?

```r
temp <- cbind(est,sample_data(M.f)[,"Treatment"])
```

```
## Error in `[.data.frame`(data.frame(x), i, j, drop = FALSE): undefined columns selected
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
## <environment: 0x10619c4a0>
## Methods may be defined for arguments: x
## Use  showMethods("t")  for currently available ones.
```

```r
dunn.test(temp[,1],temp[,2])
```

```
## Error in dunn.test(temp[, 1], temp[, 2]): object 'temp' not found
```

**Beta diversity using NMDS with Bray-Curtis as distance metric**

```r
set.seed(2)
GP.ord.BC <- ordinate(M.f, "NMDS", "bray", k=2, trymax=100)
```

```
## Square root transformation
## Wisconsin double standardization
## Run 0 stress 0.1447798 
## Run 1 stress 0.186064 
## Run 2 stress 0.1455367 
## Run 3 stress 0.1727763 
## Run 4 stress 0.1813939 
## Run 5 stress 0.1491938 
## Run 6 stress 0.1771572 
## Run 7 stress 0.1617626 
## Run 8 stress 0.1715673 
## Run 9 stress 0.1556523 
## Run 10 stress 0.1850588 
## Run 11 stress 0.1778122 
## Run 12 stress 0.1610637 
## Run 13 stress 0.1614822 
## Run 14 stress 0.1813868 
## Run 15 stress 0.1461509 
## Run 16 stress 0.1624909 
## Run 17 stress 0.1563934 
## Run 18 stress 0.1742567 
## Run 19 stress 0.1638725 
## Run 20 stress 0.1889055 
## Run 21 stress 0.1717253 
## Run 22 stress 0.1563917 
## Run 23 stress 0.1453411 
## Run 24 stress 0.1463574 
## Run 25 stress 0.1467163 
## Run 26 stress 0.1890107 
## Run 27 stress 0.1744679 
## Run 28 stress 0.1460941 
## Run 29 stress 0.1692203 
## Run 30 stress 0.1819589 
## Run 31 stress 0.1809971 
## Run 32 stress 0.1512717 
## Run 33 stress 0.1501406 
## Run 34 stress 0.1487589 
## Run 35 stress 0.1619449 
## Run 36 stress 0.174661 
## Run 37 stress 0.1858087 
## Run 38 stress 0.1461617 
## Run 39 stress 0.1766594 
## Run 40 stress 0.1524862 
## Run 41 stress 0.1658065 
## Run 42 stress 0.1745777 
## Run 43 stress 0.1652486 
## Run 44 stress 0.146165 
## Run 45 stress 0.4086956 
## Run 46 stress 0.1459195 
## Run 47 stress 0.1491701 
## Run 48 stress 0.1687325 
## Run 49 stress 0.151588 
## Run 50 stress 0.1772449 
## Run 51 stress 0.175547 
## Run 52 stress 0.173052 
## Run 53 stress 0.1745394 
## Run 54 stress 0.1814998 
## Run 55 stress 0.1761026 
## Run 56 stress 0.1802628 
## Run 57 stress 0.1749053 
## Run 58 stress 0.1770204 
## Run 59 stress 0.1656654 
## Run 60 stress 0.1654449 
## Run 61 stress 0.1574212 
## Run 62 stress 0.1670429 
## Run 63 stress 0.1552649 
## Run 64 stress 0.1736454 
## Run 65 stress 0.1657095 
## Run 66 stress 0.1855619 
## Run 67 stress 0.1683423 
## Run 68 stress 0.1790856 
## Run 69 stress 0.1677458 
## Run 70 stress 0.1746703 
## Run 71 stress 0.1740228 
## Run 72 stress 0.1554371 
## Run 73 stress 0.1751445 
## Run 74 stress 0.1862341 
## Run 75 stress 0.1455371 
## Run 76 stress 0.1563833 
## Run 77 stress 0.1681607 
## Run 78 stress 0.1696041 
## Run 79 stress 0.1792227 
## Run 80 stress 0.1478134 
## Run 81 stress 0.1788593 
## Run 82 stress 0.1478134 
## Run 83 stress 0.1739068 
## Run 84 stress 0.1832753 
## Run 85 stress 0.1747838 
## Run 86 stress 0.164247 
## Run 87 stress 0.1781949 
## Run 88 stress 0.1672539 
## Run 89 stress 0.1800106 
## Run 90 stress 0.1809296 
## Run 91 stress 0.1753395 
## Run 92 stress 0.1501244 
## Run 93 stress 0.1667888 
## Run 94 stress 0.1809751 
## Run 95 stress 0.1769462 
## Run 96 stress 0.1603592 
## Run 97 stress 0.1556344 
## Run 98 stress 0.1871302 
## Run 99 stress 0.1707339 
## Run 100 stress 0.1490922 
## *** No convergence -- monoMDS stopping criteria:
##      1: no. of iterations >= maxit
##     98: stress ratio > sratmax
##      1: scale factor of the gradient < sfgrmin
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
## Stress:     0.1447798 
## Stress type 1, weak ties
## No convergent solutions - best solution after 100 tries
## Scaling: centring, PC rotation, halfchange scaling 
## Species: expanded scores based on 'wisconsin(sqrt(veganifyOTU(physeq)))'
```

```r
color = c("Treatment")
shape = c("Dog")
title=c("NMDS of 16S microbiome,Bray-Curtis distance,k=2")
MDS = plot_ordination(M.f, GP.ord.BC, color = color,shape=shape, 
		title = title)
```

```
## Warning in plot_ordination(M.f, GP.ord.BC, color = color, shape = shape, :
## Color variable was not found in the available data you provided.No color
## mapped.
```

```
## Warning in plot_ordination(M.f, GP.ord.BC, color = color, shape = shape, :
## Shape variable was not found in the available data you provided.No shape
## mapped.
```

```r
MDS.1  = MDS +theme(axis.text=element_text(size=16, face="bold"),
				axis.title=element_text(size=18,face="bold"), legend.title=element_text(size=14))+
		theme_bw()+labs(color=color, shape=shape)+geom_point(size=5)

pdf(paste0(outDir,"/NMDS_Dogs_treatment_Bray_Curtis.pdf"),8,5)
```

```
## Error in pdf(paste0(outDir, "/NMDS_Dogs_treatment_Bray_Curtis.pdf"), 8, : cannot open file '/Users/katielennard/Documents/Academic/Postdoc/Projects/cbio_16S_pipeline_certification/github_repo/R_downstream/NMDS_Dogs_treatment_Bray_Curtis.pdf'
```

```r
MDS.1
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15-1.png)

```r
dev.off()
```

```
## RStudioGD 
##         2
```
**Beta diversity using NMDS with Unifrac as distance metric**

```r
GP.ord.U <- ordinate(M.f, "NMDS", "unifrac")
```

```
## Warning in UniFrac(physeq, ...): Randomly assigning root as -- OTU_51 -- in
## the phylogenetic tree in the data you provided.
```

```
## Run 0 stress 0.1533566 
## Run 1 stress 0.1780643 
## Run 2 stress 0.1756972 
## Run 3 stress 0.1530112 
## ... New best solution
## ... Procrustes: rmse 0.01461334  max resid 0.1036303 
## Run 4 stress 0.1737196 
## Run 5 stress 0.1838161 
## Run 6 stress 0.1590199 
## Run 7 stress 0.171394 
## Run 8 stress 0.1710265 
## Run 9 stress 0.187316 
## Run 10 stress 0.172126 
## Run 11 stress 0.1532896 
## ... Procrustes: rmse 0.0151847  max resid 0.1036364 
## Run 12 stress 0.1629224 
## Run 13 stress 0.173194 
## Run 14 stress 0.1859363 
## Run 15 stress 0.1863661 
## Run 16 stress 0.1667607 
## Run 17 stress 0.1813195 
## Run 18 stress 0.1691737 
## Run 19 stress 0.173933 
## Run 20 stress 0.1792837 
## *** No convergence -- monoMDS stopping criteria:
##     20: stress ratio > sratmax
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
## Stress:     0.1530112 
## Stress type 1, weak ties
## No convergent solutions - best solution after 20 tries
## Scaling: centring, PC rotation 
## Species: scores missing
```

```r
color = c("Treatment")
shape = c("Dog")

title=c("NMDS of 16S microbiome, Unifrac distance, k=2")

MDS = plot_ordination(M.f, GP.ord.U, color = color, shape=shape, 
		title = title)
```

```
## Warning in plot_ordination(M.f, GP.ord.U, color = color, shape = shape, :
## Color variable was not found in the available data you provided.No color
## mapped.
```

```
## Warning in plot_ordination(M.f, GP.ord.U, color = color, shape = shape, :
## Shape variable was not found in the available data you provided.No shape
## mapped.
```

```r
MDS.1  = MDS +theme(axis.text=element_text(size=16, face="bold"),
				axis.title=element_text(size=18,face="bold"), legend.title=element_text(size=14))+
		theme_bw()+labs(color=color)+geom_point(size=5)

pdf(paste0(outDir,"/NMDS_Dogs_treatment_Unifrac.pdf"),8,5)
```

```
## Error in pdf(paste0(outDir, "/NMDS_Dogs_treatment_Unifrac.pdf"), 8, 5): cannot open file '/Users/katielennard/Documents/Academic/Postdoc/Projects/cbio_16S_pipeline_certification/github_repo/R_downstream/NMDS_Dogs_treatment_Unifrac.pdf'
```

```r
MDS.1
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16-1.png)

```r
dev.off()
```

```
## RStudioGD 
##         2
```
**Create a heatmap of taxa merged at the lowest available taxonomic level**

```r
M.phy <- tax_glom.kv(M.f)#this function is available in the 'microbiome_custom_functions.R' script loaded at the beginning of this script
```

```
## [1] "Removing phylogenetic tree"
## [1] "There are now 143 merged taxa"
```

```r
ntaxa(M.phy)
```

```
## [1] 143
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

```
## Error in pdf(filename): cannot open file '/Users/katielennard/Documents/Academic/Postdoc/Projects/cbio_16S_pipeline_certification/github_repo/R_downstream/cbio_cert_heatmap_merged_taxa.pdf'
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
## Error in `[.data.frame`(data.frame(x), i, j, drop = FALSE): undefined columns selected
```

```r
print(barplot)
```

```
## function (height, ...) 
## UseMethod("barplot")
## <bytecode: 0x10481cb80>
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
## Error in ntaxa(Mraw.f): object 'Mraw.f' not found
```

```r
MGS=make_metagenomeSeq(Mraw.f)
```

```
## Error in taxa_are_rows(physeq): object 'Mraw.f' not found
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
## Error in `[.data.frame`(data.frame(x), i, j, drop = FALSE): undefined columns selected
```

```r
phy.temp <- prune_samples(sub.index, M.f)
```

```
## Error in (function (classes, fdef, mtable) : unable to find an inherited method for function 'prune_samples' for signature '"integer", "phyloseq"'
```

```r
nsamples(phy.temp)
```

```
## Error in nsamples(phy.temp): object 'phy.temp' not found
```

```r
RF.k(data = phy.temp, var = "Dog", ntree=10000, cv.fold=10, outDir = outDir, Nfeatures.validation = 3)
```

```
## Warning in file(file, ifelse(append, "a", "w")): cannot open
## file '/Users/katielennard/Documents/Academic/Postdoc/Projects/
## cbio_16S_pipeline_certification/github_repo/R_downstream/
## RFDog_results_10000_10_3_212.txt': No such file or directory
```

```
## Error in file(file, ifelse(append, "a", "w")): cannot open the connection
```

The class error rates are 0% (even one OTU enough to discriminate between Dog G and B?)

What if we used merged OTUs?

```r
merged.phy <- tax_glom.kv(phy.temp)
```

```
## Error in otu_table(physeq): object 'phy.temp' not found
```

```r
RF.k(data = merged.phy, var = "Dog", ntree=10000, cv.fold=10, outDir = outDir, Nfeatures.validation = 3, descriptor = "merged_OTUs") #for details on RF.k() see microbiome_custom_functions.R file
```

```
## Error in file(file, ifelse(append, "a", "w")): cannot open the connection
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
## Error in sample_data(physeq): object 'merged.phy' not found
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
## Error in sample_data(physeq): object 'phy.temp' not found
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
## R version 3.3.3 (2017-03-06)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## Running under: macOS Sierra 10.12.6
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] knitr_1.17           gridExtra_2.3        dplyr_0.7.4         
##  [4] metagenomeSeq_1.16.0 RColorBrewer_1.1-2   glmnet_2.0-13       
##  [7] foreach_1.4.3        Matrix_1.2-11        limma_3.30.13       
## [10] matrixStats_0.52.2   psych_1.7.8          corrplot_0.77       
## [13] dunn.test_1.3.4      ROCR_1.0-7           gplots_3.0.1        
## [16] randomForest_4.6-12  fifer_1.1            MASS_7.3-47         
## [19] Rmisc_1.5            plyr_1.8.4           ggplot2_2.2.1       
## [22] vegan_2.4-4          lattice_0.20-35      permute_0.9-4       
## [25] NMF_0.20.6           Biobase_2.34.0       BiocGenerics_0.20.0 
## [28] cluster_2.0.6        rngtools_1.2.4       pkgmaker_0.22       
## [31] registry_0.3         phyloseq_1.19.1     
## 
## loaded via a namespace (and not attached):
##  [1] nlme_3.1-131          bitops_1.0-6          doParallel_1.0.11    
##  [4] tools_3.3.3           backports_1.1.1       R6_2.2.2             
##  [7] rpart_4.1-11          KernSmooth_2.23-15    Hmisc_4.0-3          
## [10] lazyeval_0.2.0        mgcv_1.8-22           colorspace_1.3-2     
## [13] ade4_1.7-8            nnet_7.3-12           mnormt_1.5-5         
## [16] compiler_3.3.3        htmlTable_1.9         sandwich_2.4-0       
## [19] labeling_0.3          caTools_1.17.1        scales_0.5.0         
## [22] checkmate_1.8.4       mvtnorm_1.0-6         stringr_1.2.0        
## [25] digest_0.6.12         foreign_0.8-69        XVector_0.14.1       
## [28] base64enc_0.1-3       pkgconfig_2.0.1       htmltools_0.3.6      
## [31] plotrix_3.6-6         highr_0.6             maps_3.2.0           
## [34] htmlwidgets_0.9       rlang_0.1.2           bindr_0.1            
## [37] zoo_1.8-0             jsonlite_1.5          gtools_3.5.0         
## [40] acepack_1.4.1         magrittr_1.5          modeltools_0.2-21    
## [43] Formula_1.2-2         dotCall64_0.9-04      biomformat_1.2.0     
## [46] Rcpp_0.12.13          munsell_0.4.3         S4Vectors_0.12.2     
## [49] ape_4.1               stringi_1.1.5         multcomp_1.4-7       
## [52] zlibbioc_1.20.0       rhdf5_2.18.0          grid_3.3.3           
## [55] strucchange_1.5-1     gdata_2.18.0          Biostrings_2.42.1    
## [58] splines_3.3.3         multtest_2.30.0       igraph_1.1.2         
## [61] party_1.2-3           reshape2_1.4.2        codetools_0.2-15     
## [64] stats4_3.3.3          glue_1.1.1            evaluate_0.10.1      
## [67] latticeExtra_0.6-28   data.table_1.10.4     spam_2.1-1           
## [70] gtable_0.2.0          assertthat_0.2.0      gridBase_0.4-7       
## [73] coin_1.2-1            xtable_1.8-2          survival_2.41-3      
## [76] randomForestSRC_2.5.0 tibble_1.3.4          iterators_1.0.8      
## [79] IRanges_2.8.2         bindrcpp_0.2          fields_9.0           
## [82] TH.data_1.0-8
```
