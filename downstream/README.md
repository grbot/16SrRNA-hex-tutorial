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

Import data and create phyloseq object
--------------------------------------

**Import BIOM file (generated in QIIME) into a phyloseq object**

    library(phyloseq)
    library(ggplot2)
    library(gridExtra)
    library(dunn.test)

**Import custom functions used in script**

    source("/home/gerrit/workspace/amw/src/microbiome_custom_functions.R")

    ## Loading required package: pkgmaker

    ## Loading required package: registry

    ## 
    ## Attaching package: 'pkgmaker'

    ## The following object is masked from 'package:base':
    ## 
    ##     isNamespaceLoaded

    ## Loading required package: rngtools

    ## Loading required package: cluster

    ## NMF - BioConductor layer [OK] | Shared memory capabilities [NO: bigmemory] | Cores 3/4

    ##   To enable shared memory capabilities, try: install.extras('
    ## NMF
    ## ')

    ## Loading required package: permute

    ## Loading required package: lattice

    ## This is vegan 2.4-2

    ## 
    ## Attaching package: 'psych'

    ## The following objects are masked from 'package:ggplot2':
    ## 
    ##     %+%, alpha

    ## matrixStats v0.51.0 (2016-10-08) successfully loaded. See ?matrixStats for help.

    ## 
    ## Attaching package: 'matrixStats'

    ## The following objects are masked from 'package:Biobase':
    ## 
    ##     anyMissing, rowMedians

    ## Loading required package: MASS

    ## 
    ## Attaching package: 'fifer'

    ## The following object is masked from 'package:Biobase':
    ## 
    ##     contents

    ## Loading required package: limma

    ## 
    ## Attaching package: 'limma'

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     plotMA

    ## Loading required package: glmnet

    ## Loading required package: Matrix

    ## Loading required package: foreach

    ## Loaded glmnet 2.0-5

    ## Loading required package: RColorBrewer

    ## randomForest 4.6-12

    ## Type rfNews() to see new features/changes/bug fixes.

    ## 
    ## Attaching package: 'randomForest'

    ## The following object is masked from 'package:psych':
    ## 
    ##     outlier

    ## The following object is masked from 'package:Biobase':
    ## 
    ##     combine

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     combine

    ## The following object is masked from 'package:gridExtra':
    ## 
    ##     combine

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     margin

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:randomForest':
    ## 
    ##     combine

    ## The following object is masked from 'package:MASS':
    ## 
    ##     select

    ## The following object is masked from 'package:matrixStats':
    ## 
    ##     count

    ## The following object is masked from 'package:Biobase':
    ## 
    ##     combine

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following object is masked from 'package:gridExtra':
    ## 
    ##     combine

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    ## Loading required package: gplots

    ## 
    ## Attaching package: 'gplots'

    ## The following object is masked from 'package:stats':
    ## 
    ##     lowess

**Set the working directory and import data**

    setwd("/home/gerrit/scratch/amw")
    inDir <- getwd()
    outDir <- paste0(getwd(),"/downstream") # Specify output directory
    phy <- import_biom(BIOMfilename = paste0(inDir,"/otus_table.tax.biom"), 
            verbose = TRUE)#
    ntaxa(phy)#(number of OTUs)

    ## [1] 181

    sample_names(phy) <- sub("\\/1","",sample_names(phy))#remove "/1" from filenames
    # Add phylogenetic tree (.tre file generated in QIIME)
    tree <- read_tree_greengenes(paste0(inDir,"/otus_repsetOUT_aligned_pfiltered.tre"))
    # Merge phy and tree
    phy <- merge_phyloseq(phy,tree)

**Data cleanup**

    colnames(tax_table(phy))

    ## [1] "Rank1" "Rank2" "Rank3" "Rank4" "Rank5" "Rank6" "Rank7"

    colnames(tax_table(phy)) <-  c("Kingdom", "Phylum" , "Class" , "Order" , "Family" , "Genus", "Species")# e.g. Replace "Rank1" with "Kingdom"
    # Clean taxonomic annotations, at the moment they are for example 'k__Bacteria'; 'p_Firmicutes' - remove k__ and p__ ...
    tax_table(phy)[,"Kingdom"] <- sub("k__","",tax_table(phy)[,"Kingdom"])
    tax_table(phy)[,"Phylum"] <- sub("p__","",tax_table(phy)[,"Phylum"])
    tax_table(phy)[,"Class"] <- sub("c__","",tax_table(phy)[,"Class"])
    tax_table(phy)[,"Order"] <- sub("o__","",tax_table(phy)[,"Order"])
    tax_table(phy)[,"Family"] <- sub("f__","",tax_table(phy)[,"Family"])
    tax_table(phy)[,"Genus"] <- sub("g__","",tax_table(phy)[,"Genus"])
    tax_table(phy)[,"Species"] <- sub("s__","",tax_table(phy)[,"Species"])

**Need to filter out unclassified OTUs otherwise custom functions will
fail**

    t= which(is.na(tax_table(phy)[,"Phylum"]))
    tax_table(phy) = tax_table(phy)[-t,]

**Import metadata and merge with phyloseq object**

    meta <-  read.table(paste0(inDir,"/practice.dataset1.metadata.tsv"), sep = "\t", header =TRUE, row.names=1)
    head(meta)

    ##       Dog Treatment
    ## Dog1    B         2
    ## Dog2    G         3
    ## Dog3    K         3
    ## Dog8    B         4
    ## Dog9    G         0
    ## Dog10   K         4

    rownames(meta)

    ##  [1] "Dog1"  "Dog2"  "Dog3"  "Dog8"  "Dog9"  "Dog10" "Dog15" "Dog16"
    ##  [9] "Dog17" "Dog22" "Dog23" "Dog24" "Dog29" "Dog30" "Dog31"

    head(sample_names(phy))

    ## [1] "Dog10" "Dog15" "Dog16" "Dog17" "Dog1"  "Dog22"

    length(sample_names(phy))#15

    ## [1] 15

    length(rownames(meta))#15 (check if same number of samples in .biom file and metadatafile)

    ## [1] 15

    length(intersect(rownames(meta),sample_names(phy)))#15 (check that the sample names match in all cases)

    ## [1] 15

    sample_data(phy) <- meta # Assign the metadata to the phyloseq object 'phy' (phyloseq will put these in the right order)
    nsamples(phy)

    ## [1] 15

    str(sample_data(phy)) # Need to change treatment column to factor variable

    ## 'data.frame':    15 obs. of  2 variables:
    ## Formal class 'sample_data' [package "phyloseq"] with 4 slots
    ##   ..@ .Data    :List of 2
    ##   .. ..$ : Factor w/ 3 levels "B","G","K": 3 1 2 3 1 1 2 3 1 2 ...
    ##   .. ..$ : int  4 1 4 0 2 3 1 2 0 3 ...
    ##   ..@ names    : chr  "Dog" "Treatment"
    ##   ..@ row.names: chr  "Dog10" "Dog15" "Dog16" "Dog17" ...
    ##   ..@ .S3Class : chr "data.frame"

    sample_data(phy)[,"Treatment"] <- as.numeric(unlist(sample_data(phy)[,"Treatment"]))

**Save phyloseq object as an .RData file**

    save(phy, file = paste0(outDir,"/dog_stool.RData")) # Save annotated object as a .RData object
    load(paste0(outDir,"/dog_stool.RData"))

Explore number of reads per sample, make rarefaction curves and filter data as necessary
----------------------------------------------------------------------------------------

**Explore number of reads per sample**

    reads <- sample_sums(phy)
    length(which(reads<5000))

    ## [1] 0

    raremax <- min(reads)
    raremax

    ## [1] 63966

    rarecurve(t(otu_table(phy)), step = 100, sample = raremax,xlab = "number of reads/sample", ylab = "number of OTUs",
            label = FALSE, xlim = c(0,100000))

![](README_files/figure-markdown_strict/unnamed-chunk-8-1.png)

All samples have sufficient sequencing depth for inclusion in downstream
analyses. The vertical line in the above plot indicates the sample with
the lowest number of reads. Now we will scale data to account for
differences in the number of reads/sample and filter rare OTUs that are
not of biological interest for the purpose of this analysis (e.g. occurs
only in one sample). **Standardize abundances to median sequence depth**

    total = median(sample_sums(phy))
    standf = function(x, t=total) round(t * (x / sum(x)))
    M.std = transform_sample_counts(phy, standf)

**Apply mild OTU filter**

Select OTUs where the rowsum for that OTU has at least 20% of samples
with a count of 10 each OR where that OTU &gt; 0.001% of the total
median count (for cases where the minority of samples may have high
counts of a rare OTU)

    M.f = filter_taxa(M.std,function(x) sum(x > 10) > (0.02*length(x)) | sum(x) > 0.001*total, TRUE)
    ntaxa(M.f)

    ## [1] 138

**Basic exploratory plots: alpha- and beta-diversity, barplots, heatmap**
-------------------------------------------------------------------------

**Alpha diversity by dog**

    p <- plot_richness(M.std,x = "Dog",color = "Treatment",measures=c("Shannon"), 
            title = paste0("Standardized to total reads, N=",nsamples(M.std)))+theme(axis.text=element_text(size=16, face="bold"),
                    axis.title=element_text(size=16,face="bold"))+geom_point(size=5)
    p

![](README_files/figure-markdown_strict/unnamed-chunk-11-1.png)

    pdf(paste0(outDir,"/alpha_diversity_by_dog_treatment.pdf"))
    p
    dev.off()

    ## PNG 
    ##   2

Is there a significant difference in alpha diversity between dogs
irrespective of treatment?

    est <- estimate_richness(M.f, split = TRUE, measures = c("Shannon"))
    temp <- cbind(est,sample_data(M.f)[,"Dog"])
    head(temp)

    ##        Shannon Dog
    ## Dog10 2.893430   K
    ## Dog15 2.566590   B
    ## Dog16 2.930422   G
    ## Dog17 2.850480   K
    ## Dog1  3.216201   B
    ## Dog22 2.538645   B

    t <- kruskal.test(temp[,1]~temp[,2])
    t

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  temp[, 1] by temp[, 2]
    ## Kruskal-Wallis chi-squared = 4.94, df = 2, p-value = 0.08458

    dunn.test(temp[,1],temp[,2])#post-hoc testing to see which dogs are different

    ##   Kruskal-Wallis rank sum test
    ## 
    ## data: x and group
    ## Kruskal-Wallis chi-squared = 4.94, df = 2, p-value = 0.08
    ## 
    ## 
    ##                            Comparison of x by group                            
    ##                                 (No adjustment)                                
    ## Col Mean-|
    ## Row Mean |          B          G
    ## ---------+----------------------
    ##        G |  -1.767766
    ##          |     0.0385
    ##          |
    ##        K |   0.282842   2.050609
    ##          |     0.3886     0.0202

Dog G has significantly higher alpha diversity than dogs K and B
irrespective of treatment

**Alpha diversity by treatment**

    p <- plot_richness(M.std,x = "Treatment",color = "Dog",measures=c("Shannon"), 
                    title = paste0("Standardized to total reads, N=",nsamples(M.std)))+theme(axis.text=element_text(size=16, face="bold"),
                    axis.title=element_text(size=16,face="bold"))+geom_point(size=5)
    p

![](README_files/figure-markdown_strict/unnamed-chunk-13-1.png)

    pdf(paste0(outDir,"/alpha_diversity_by_treatment_dog.pdf"))
    p
    dev.off()

    ## PNG 
    ##   2

Are there significant differences in alpha diversity by treatment?

    temp <- cbind(est,sample_data(M.f)[,"Treatment"])
    head(temp)

    ##        Shannon Treatment
    ## Dog10 2.893430         4
    ## Dog15 2.566590         1
    ## Dog16 2.930422         4
    ## Dog17 2.850480         0
    ## Dog1  3.216201         2
    ## Dog22 2.538645         3

    t <- kruskal.test(temp[,1]~temp[,2])
    t

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  temp[, 1] by temp[, 2]
    ## Kruskal-Wallis chi-squared = 1.6, df = 4, p-value = 0.8088

    dunn.test(temp[,1],temp[,2])

    ##   Kruskal-Wallis rank sum test
    ## 
    ## data: x and group
    ## Kruskal-Wallis chi-squared = 1.6, df = 4, p-value = 0.81
    ## 
    ## 
    ##                            Comparison of x by group                            
    ##                                 (No adjustment)                                
    ## Col Mean-|
    ## Row Mean |          4          1          0          2
    ## ---------+--------------------------------------------
    ##        1 |   0.273861
    ##          |     0.3921
    ##          |
    ##        0 |  -0.730296  -1.004158
    ##          |     0.2326     0.1577
    ##          |
    ##        2 |  -0.182574  -0.456435   0.547722
    ##          |     0.4276     0.3240     0.2919
    ##          |
    ##        3 |  -0.730296  -1.004158   0.000000  -0.547722
    ##          |     0.2326     0.1577     0.5000     0.2919

**Beta diversity using NMDS with Bray-Curtis as distance metric**

    set.seed(2)
    GP.ord.BC <- ordinate(M.f, "NMDS", "bray", k=2, trymax=100) # stress=0.09

    ## Square root transformation
    ## Wisconsin double standardization
    ## Run 0 stress 0.07130578 
    ## Run 1 stress 0.0709321 
    ## ... New best solution
    ## ... Procrustes: rmse 0.009855629  max resid 0.0279429 
    ## Run 2 stress 0.0721892 
    ## Run 3 stress 0.1100664 
    ## Run 4 stress 0.07329299 
    ## Run 5 stress 0.3620402 
    ## Run 6 stress 0.07218952 
    ## Run 7 stress 0.07398944 
    ## Run 8 stress 0.2138891 
    ## Run 9 stress 0.07241017 
    ## Run 10 stress 0.07093206 
    ## ... New best solution
    ## ... Procrustes: rmse 5.022355e-05  max resid 0.0001204355 
    ## ... Similar to previous best
    ## Run 11 stress 0.1190393 
    ## Run 12 stress 0.07093231 
    ## ... Procrustes: rmse 0.0003115158  max resid 0.0007574303 
    ## ... Similar to previous best
    ## Run 13 stress 0.07218918 
    ## Run 14 stress 0.1100655 
    ## Run 15 stress 0.07093235 
    ## ... Procrustes: rmse 0.0003326004  max resid 0.0008081068 
    ## ... Similar to previous best
    ## Run 16 stress 0.1100645 
    ## Run 17 stress 0.07241209 
    ## Run 18 stress 0.07218928 
    ## Run 19 stress 0.07399004 
    ## Run 20 stress 0.07241308 
    ## *** Solution reached

    color = c("Treatment")
    shape = c("Dog")
    title=c("NMDS of 16S microbiome,Bray-Curtis distance,k=2")
    MDS = plot_ordination(M.f, GP.ord.BC, color = color,shape=shape, 
            title = title)
    MDS.1  = MDS +theme(axis.text=element_text(size=16, face="bold"),
                    axis.title=element_text(size=18,face="bold"), legend.title=element_text(size=14))+
            theme_bw()+labs(color=color, shape=shape)+geom_point(size=5)

    MDS.1

![](README_files/figure-markdown_strict/unnamed-chunk-15-1.png)

    pdf(paste0(outDir,"/NMDS_Dogs_tretment_Bray_Curtis.pdf"),8,5)
    MDS.1
    dev.off()

    ## PNG 
    ##   2

**Beta diversity using NMDS with Unifrac as distance metric**

    GP.ord.U <- ordinate(M.f, "NMDS", "unifrac")#stress=0.08

    ## Warning in UniFrac(physeq, ...): Randomly assigning root as -- OTU_179 --
    ## in the phylogenetic tree in the data you provided.

    ## Run 0 stress 0.1020741 
    ## Run 1 stress 0.1020741 
    ## ... New best solution
    ## ... Procrustes: rmse 7.693494e-05  max resid 0.0001512306 
    ## ... Similar to previous best
    ## Run 2 stress 0.1005345 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1036481  max resid 0.2806968 
    ## Run 3 stress 0.09338723 
    ## ... New best solution
    ## ... Procrustes: rmse 0.08563969  max resid 0.2850397 
    ## Run 4 stress 0.1121452 
    ## Run 5 stress 0.09338727 
    ## ... Procrustes: rmse 8.017275e-05  max resid 0.0001472236 
    ## ... Similar to previous best
    ## Run 6 stress 0.1521734 
    ## Run 7 stress 0.09338724 
    ## ... Procrustes: rmse 3.756729e-05  max resid 7.178452e-05 
    ## ... Similar to previous best
    ## Run 8 stress 0.1005342 
    ## Run 9 stress 0.1343305 
    ## Run 10 stress 0.1005345 
    ## Run 11 stress 0.1020741 
    ## Run 12 stress 0.09338727 
    ## ... Procrustes: rmse 8.519511e-05  max resid 0.00015642 
    ## ... Similar to previous best
    ## Run 13 stress 0.09338723 
    ## ... New best solution
    ## ... Procrustes: rmse 1.734821e-05  max resid 3.107499e-05 
    ## ... Similar to previous best
    ## Run 14 stress 0.1020741 
    ## Run 15 stress 0.1005343 
    ## Run 16 stress 0.1121452 
    ## Run 17 stress 0.1005343 
    ## Run 18 stress 0.09338724 
    ## ... Procrustes: rmse 1.616961e-05  max resid 3.354203e-05 
    ## ... Similar to previous best
    ## Run 19 stress 0.09338733 
    ## ... Procrustes: rmse 5.409304e-05  max resid 0.0001203802 
    ## ... Similar to previous best
    ## Run 20 stress 0.112145 
    ## *** Solution reached

    color = c("Treatment")
    shape = c("Dog")

    title=c("NMDS of 16S microbiome, Unifrac distance,k=2")

    MDS = plot_ordination(M.f, GP.ord.U, color = color, shape=shape, 
            title = title)
    MDS.1  = MDS +theme(axis.text=element_text(size=16, face="bold"),
                    axis.title=element_text(size=18,face="bold"), legend.title=element_text(size=14))+
            theme_bw()+labs(color=color)+geom_point(size=5)
    MDS.1

![](README_files/figure-markdown_strict/unnamed-chunk-16-1.png)

    pdf(paste0(outDir,"/NMDS_Dogs_treatment_Bray_Curtis.pdf"),8,5)
    MDS.1
    dev.off()

    ## PNG 
    ##   2

**Create a heatmap of taxa merged at the lowest available taxonomic
level**

    M.phy <- tax_glom.kv(M.f) # This function is available in the 'microbiome_custom_functions.R' script loaded at the beginning of this script

    ## [1] "Removing phylogenetic tree"
    ## [1] "There are now 58 merged taxa"

    ntaxa(M.phy)

    ## [1] 58

    filename <- c("heatmap_merged_taxa")
    main <- paste("Merged taxa, Bray-Curtis distance")
    f = paste0(outDir,filename,".pdf")
    # Color specification for column annotations above heatmap:
    D.cols = c("B"="#CC79A7","G"="#56B4E9","K"="#F0E442")
    colours = list(Dog=D.cols)

    # Create distance matrix and calculate tree:
    set.seed(2)
    diss <- distance(M.phy,method = "bray", type = "samples")
    clust.res<-hclust(diss)
    sample.order = clust.res$order
    # Heatmap is output to file (the heatmap.k function can be found in the 'microbiome_custom_functions.R' script)
    hm = heatmap.k(physeq= M.phy, annot.cols = c(1,2), main = main,filename = f,colours=colours,Colv = sample.order,labrow = TRUE, cexCol = 2)  

    ## [1] "including all samples"
    ## [1] "including all otus"

![](README_files/figure-markdown_strict/unnamed-chunk-17-1.png)

    print(hm)

    ## $Rowv
    ## 'dendrogram' with 2 branches and 58 members total, at height 31.61713 
    ## 
    ## $rowInd
    ##  [1] 16 44 30 50 12 53 39 42 28 33 35 18 21 19 31  7 20 38 47  3 36  9 43
    ## [24] 27  6 14 55  4 29 54 49 48 56  5  8  2 22 17 26 41 37 24 34 57 13 58
    ## [47] 23 11  1 32 46 51 15 45 52 40 25 10
    ## 
    ## $colInd
    ##  [1] 13  9  2  6  5 10 15  3  7 11  1 14 12  4  8

**Barplots by dog**
-------------------

    level = "Genus"
    count = 500
    perc = 0.25
    # Barplot will be written to file (the bar.plots function can be found in the 'microbiome_custom_functions.R' script)
    barplot = bar.plots(physeq = M.std,cat = "Dog",level = level, count = count, perc = perc, outDir=outDir, 
            filen = 'Barplots_by_Dog')
    print(barplot)

![](README_files/figure-markdown_strict/unnamed-chunk-18-1.png)

Detect taxa/OTUs that differ significantly by Dog
-------------------------------------------------

convert phyloseq object to metagenomeSeq obj. NB use raw data not
standardized:

    Mraw.f = filter_taxa(phy,function(x) sum(x > 10) > (0.02*length(x)) | sum(x) > 0.001*total, TRUE)
    ntaxa(Mraw.f)

    ## [1] 138

    MGS=make_metagenomeSeq(Mraw.f)

    ## Default value being used.

    MGS

    ## MRexperiment (storageMode: environment)
    ## assayData: 138 features, 15 samples 
    ##   element names: counts 
    ## protocolData: none
    ## phenoData
    ##   sampleNames: Dog10 Dog15 ... Dog9 (15 total)
    ##   varLabels: Dog Treatment
    ##   varMetadata: labelDescription
    ## featureData
    ##   featureNames: OTU_72 OTU_77 ... OTU_47 (138 total)
    ##   fvarLabels: OTUname
    ##   fvarMetadata: labelDescription
    ## experimentData: use 'experimentData(object)'
    ## Annotation:

**Use Random forests analysis to detect taxa that are good predictors of
Dog**

Example used: Dog G vs. Dog B (all treatment points)

    sub.index <- sample_names(M.f)[sample_data(M.f)[,"Dog"] != "K"]
    phy.temp <- prune_samples(sub.index, M.f)
    nsamples(phy.temp)

    ## [1] 10

    library(randomForest)
    library(dplyr)
    RF.k(data = phy.temp, var = "Dog", ntree=10000, cv.fold=10, outDir = outDir, Nfeatures.validation = 3)

    ## [1] "0 samples did not have response variable data, removing these..."
    ## [1] "Data set size:  10 samples with 5 and 5 samples per class"
    ## [1] "Cross-validated error rates associated with stepwise reduction of features:"
    ## 138  69  34  17   9   4   2   1 
    ## 0.0 0.0 0.0 0.0 0.0 0.0 0.2 0.2

    ## [1] "*****************************"
    ## [1] "THE TOP 20 MOST IMPORTANT FEATURES WERE:"
    ##         predictors        B        G MeanDecreaseAccuracy MeanDecreaseGini
    ## OTU_141    OTU_141 8.853276 8.938303             9.195140       0.08360000
    ## OTU_125    OTU_125 9.852233 9.753600            10.139613       0.08078000
    ## OTU_13      OTU_13 8.654537 8.670092             8.968530       0.07870000
    ## OTU_85      OTU_85 9.234499 9.270183             9.532235       0.07676000
    ## OTU_26      OTU_26 7.610852 7.415002             7.867611       0.07600667
    ## OTU_28      OTU_28 9.269239 9.179335             9.553787       0.07532000
    ## OTU_45      OTU_45 9.163334 9.093742             9.543887       0.07484000
    ## OTU_82      OTU_82 9.697815 9.603315             9.918552       0.07482000
    ## OTU_44      OTU_44 8.965029 8.997229             9.289003       0.07447000
    ## OTU_92      OTU_92 8.255487 8.218606             8.530282       0.07438000
    ## OTU_9        OTU_9 9.519245 9.459887             9.890981       0.07430000
    ## OTU_11      OTU_11 8.910690 8.823479             9.208113       0.07396000
    ## OTU_78      OTU_78 9.283935 9.161688             9.611126       0.07384000
    ## OTU_56      OTU_56 9.681690 9.679827            10.023366       0.07378000
    ## OTU_101    OTU_101 8.796644 8.698407             9.041628       0.07350000
    ## OTU_19      OTU_19 8.805800 8.729367             9.031053       0.07338667
    ## OTU_24      OTU_24 8.008401 7.988265             8.385150       0.07324000
    ## OTU_34      OTU_34 8.668624 8.571000             8.930318       0.07220000
    ## OTU_47      OTU_47 7.273510 7.348116             7.651533       0.07202667
    ## OTU_110    OTU_110 7.982483 7.935700             8.204203       0.07160000
    ##                          tax
    ## OTU_141        Fusobacterium
    ## OTU_125         Enterococcus
    ## OTU_13         F.prausnitzii
    ## OTU_85   Erysipelotrichaceae
    ## OTU_26       Lachnospiraceae
    ## OTU_28          Turicibacter
    ## OTU_45                 Dorea
    ## OTU_82    [Mogibacteriaceae]
    ## OTU_44               Dorea.1
    ## OTU_92          Oscillospira
    ## OTU_9                P.copri
    ## OTU_11            B.producta
    ## OTU_78         Clostridiales
    ## OTU_56         Bacteroidales
    ## OTU_101       Oscillospira.1
    ## OTU_19     Lachnospiraceae.1
    ## OTU_24                 S24-7
    ## OTU_34  [Paraprevotellaceae]
    ## OTU_47     Lachnospiraceae.2
    ## OTU_110      Ruminococcaceae
    ## [1] "*****************************"

    ## [1] "Training AUC=1"
    ## [1] "Training PPV=1"
    ## [1] "Training NPV=1"

    ## [1] "*****************************"
    ## [1] "Training set classification summary if using the top 3 features only"
    ## [1] "Feature(s) selected: OTU_141" "Feature(s) selected: OTU_125"
    ## [3] "Feature(s) selected: OTU_13" 
    ##                 B        G MeanDecreaseAccuracy MeanDecreaseGini
    ## OTU_141 0.1063500 0.108525           0.09641381          1.55448
    ## OTU_13  0.1021667 0.101050           0.09129048          1.48478
    ## OTU_125 0.1358417 0.131375           0.11959905          1.47374
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

The class error rates are 0% (even one OTU enough to discriminate
between Dog G and B?)

What if we used merged OTUs?

    merged.phy <- tax_glom.kv(phy.temp)

    ## [1] "Removing phylogenetic tree"
    ## [1] "There are now 58 merged taxa"

    RF.k(data = merged.phy, var = "Dog", ntree=10000, cv.fold=10, outDir = outDir, Nfeatures.validation = 3, descriptor = "merged_OTUs")

    ## [1] "0 samples did not have response variable data, removing these..."
    ## [1] "Data set size:  10 samples with 5 and 5 samples per class"
    ## [1] "Cross-validated error rates associated with stepwise reduction of features:"
    ##  58  29  14   7   4   1 
    ## 0.0 0.0 0.0 0.0 0.0 0.1

    ## [1] "*****************************"
    ## [1] "THE TOP 20 MOST IMPORTANT FEATURES WERE:"
    ##         predictors        B        G MeanDecreaseAccuracy MeanDecreaseGini
    ## OTU_73      OTU_73 13.92293 13.90158             14.48344        0.1737167
    ## OTU_10      OTU_10 12.49247 12.47807             12.94437        0.1707314
    ## OTU_13      OTU_13 11.55367 11.68032             11.98464        0.1706067
    ## OTU_54      OTU_54 12.84943 12.81846             13.41733        0.1681267
    ## OTU_125    OTU_125 13.87000 13.90372             14.39531        0.1656083
    ## OTU_56      OTU_56 13.29844 13.39568             13.86765        0.1643100
    ## OTU_71      OTU_71 13.34526 13.26860             13.77618        0.1642600
    ## OTU_34      OTU_34 13.90986 13.99704             14.43902        0.1635067
    ## OTU_95      OTU_95 14.07960 13.96190             14.56578        0.1631314
    ## OTU_59      OTU_59 12.42655 12.48879             12.91485        0.1584114
    ## OTU_24      OTU_24 12.74488 12.67917             13.22373        0.1579200
    ## OTU_111    OTU_111 12.67367 12.64432             13.03962        0.1540067
    ## OTU_28      OTU_28 13.59523 13.63420             14.13919        0.1526267
    ## OTU_35      OTU_35 12.56198 12.62705             13.19065        0.1508400
    ## OTU_11      OTU_11 12.86726 12.84785             13.43512        0.1495933
    ## OTU_9        OTU_9 13.33599 13.39256             13.93159        0.1492267
    ## OTU_82      OTU_82 13.33169 13.25126             13.79476        0.1484067
    ## OTU_27      OTU_27 13.33507 13.27208             13.76675        0.1484000
    ## OTU_18      OTU_18 12.64827 12.65739             13.12016        0.1445700
    ## OTU_68      OTU_68 13.48077 13.49082             13.97638        0.1441900
    ##                           tax
    ## OTU_73           Oscillospira
    ## OTU_10     Enterobacteriaceae
    ## OTU_13          F.prausnitzii
    ## OTU_54    Succinivibrionaceae
    ## OTU_125          Enterococcus
    ## OTU_56          Bacteroidales
    ## OTU_71  Peptostreptococcaceae
    ## OTU_34   [Paraprevotellaceae]
    ## OTU_95                  Dorea
    ## OTU_59         [Ruminococcus]
    ## OTU_24                  S24-7
    ## OTU_111           Odoribacter
    ## OTU_28           Turicibacter
    ## OTU_35            Allobaculum
    ## OTU_11             B.producta
    ## OTU_9                 P.copri
    ## OTU_82     [Mogibacteriaceae]
    ## OTU_27           Ruminococcus
    ## OTU_18               R.gnavus
    ## OTU_68          Adlercreutzia
    ## [1] "*****************************"

    ## [1] "Training AUC=1"
    ## [1] "Training PPV=1"
    ## [1] "Training NPV=1"

    ## [1] "*****************************"
    ## [1] "Training set classification summary if using the top 3 features only"
    ## [1] "Feature(s) selected: OTU_73" "Feature(s) selected: OTU_10"
    ## [3] "Feature(s) selected: OTU_13"
    ##                 B          G MeanDecreaseAccuracy MeanDecreaseGini
    ## OTU_73 0.14227500 0.14479167           0.12855762          1.55290
    ## OTU_10 0.11312500 0.10742500           0.09887333          1.50766
    ## OTU_13 0.09850833 0.09779167           0.08807857          1.45636
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

**Differential abundance testing using MetagenomeSeq package**

Lets again compare dog G vs. dog B (merged taxa)

    colours = list(Dog=D.cols)
    a = super.fitZig.kv(physeq = merged.phy,factor = "Dog",outDir = outDir,FileName =c("1_25FC_0.2_Dog_GvsB_taxa_merged"),
    heatmap.descriptor=c("tax_annot"), main=c("Dog G vs. B, taxa merged"), subt=c("subt = FDR < ###0.05,|coeff| >= 1.25, >20%+ in either group"), ordered=TRUE, p=0.05, FC = 1.25, perc=0.2, extra.cols = c("Treatment"))

    ## [1] "0 of 10 samples were removed due to missing data"

    ## Default value being used.

    ## [1] "Dog will be modeled as a binary categorical predictor variable"

    ## Default value being used.

    ## it= 0, nll=18.47, log10(eps+1)=Inf, stillActive=58
    ## it= 1, nll=19.48, log10(eps+1)=0.01, stillActive=4
    ## it= 2, nll=19.52, log10(eps+1)=0.03, stillActive=1
    ## it= 3, nll=19.56, log10(eps+1)=0.00, stillActive=1
    ## it= 4, nll=19.64, log10(eps+1)=0.00, stillActive=0
    ## There were  25 OTUs significantly different between B vs. G that met 
    ##  threshold criteria of p 0.05 absolute FC 1.25 and percentage presence in at least one group of 20 % 
    ## [1] "writing results and model to file"
    ## [1] "/home/gerrit/scratch/amw/downstream/1_25FC_0.2_Dog_GvsB_taxa_merged_tax_annot.pdf"

![](README_files/figure-markdown_strict/unnamed-chunk-22-1.png)![](README_files/figure-markdown_strict/unnamed-chunk-22-2.png)

    ## [1] "making heatmap of results"

    print(a)

    ##         percent_positive_group0 percent_positive_group1
    ## OTU_10                      100                     100
    ## OTU_106                      80                      20
    ## OTU_18                      100                     100
    ## OTU_74                      100                      80
    ## OTU_28                      100                     100
    ## OTU_125                     100                      40
    ## OTU_71                      100                     100
    ## OTU_95                      100                     100
    ## OTU_11                      100                     100
    ## OTU_59                      100                     100
    ## OTU_70                      100                     100
    ## OTU_13                      100                     100
    ## OTU_42                       80                     100
    ## OTU_73                      100                     100
    ## OTU_27                      100                     100
    ## OTU_82                       60                     100
    ## OTU_66                       60                     100
    ## OTU_38                      100                     100
    ## OTU_9                       100                     100
    ## OTU_24                      100                     100
    ## OTU_54                       60                     100
    ## OTU_68                        0                     100
    ## OTU_64                       40                     100
    ## OTU_56                       20                     100
    ## OTU_34                       40                     100
    ##         +samples in group 0 +samples in group 1 mean_positive_group0
    ## OTU_10                    5                   5                 3645
    ## OTU_106                   4                   1                   32
    ## OTU_18                    5                   5                 1142
    ## OTU_74                    5                   4                   68
    ## OTU_28                    5                   5                 1147
    ## OTU_125                   5                   2                   10
    ## OTU_71                    5                   5                  167
    ## OTU_95                    5                   5                 1570
    ## OTU_11                    5                   5                 4034
    ## OTU_59                    5                   5                  105
    ## OTU_70                    5                   5                  119
    ## OTU_13                    5                   5                 1822
    ## OTU_42                    4                   5                   37
    ## OTU_73                    5                   5                   16
    ## OTU_27                    5                   5                  242
    ## OTU_82                    3                   5                    7
    ## OTU_66                    3                   5                    4
    ## OTU_38                    5                   5                   85
    ## OTU_9                     5                   5                 1524
    ## OTU_24                    5                   5                   20
    ## OTU_54                    3                   5                    2
    ## OTU_68                    0                   5                  NaN
    ## OTU_64                    2                   5                    3
    ## OTU_56                    1                   5                    4
    ## OTU_34                    2                   5                   10
    ##         mean_positive_group1 oddsRatio      lower       upper     fisherP
    ## OTU_10                    70    0.0000 0.00000000         Inf 1.000000000
    ## OTU_106                    1   10.9072 0.45473998 968.7617574 0.206349206
    ## OTU_18                    85    0.0000 0.00000000         Inf 1.000000000
    ## OTU_74                     6       Inf 0.02564066         Inf 1.000000000
    ## OTU_28                   241    0.0000 0.00000000         Inf 1.000000000
    ## OTU_125                    1       Inf 0.49337123         Inf 0.166666667
    ## OTU_71                    28    0.0000 0.00000000         Inf 1.000000000
    ## OTU_95                   409    0.0000 0.00000000         Inf 1.000000000
    ## OTU_11                  1279    0.0000 0.00000000         Inf 1.000000000
    ## OTU_59                    32    0.0000 0.00000000         Inf 1.000000000
    ## OTU_70                    38    0.0000 0.00000000         Inf 1.000000000
    ## OTU_13                   588    0.0000 0.00000000         Inf 1.000000000
    ## OTU_42                   201    0.0000 0.00000000  39.0005500 1.000000000
    ## OTU_73                   121    0.0000 0.00000000         Inf 1.000000000
    ## OTU_27                  1452    0.0000 0.00000000         Inf 1.000000000
    ## OTU_82                    72    0.0000 0.00000000   5.1183766 0.444444444
    ## OTU_66                    88    0.0000 0.00000000   5.1183766 0.444444444
    ## OTU_38                   598    0.0000 0.00000000         Inf 1.000000000
    ## OTU_9                  15540    0.0000 0.00000000         Inf 1.000000000
    ## OTU_24                   455    0.0000 0.00000000         Inf 1.000000000
    ## OTU_54                   104    0.0000 0.00000000   5.1183766 0.444444444
    ## OTU_68                    52    0.0000 0.00000000   0.4353226 0.007936508
    ## OTU_64                   231    0.0000 0.00000000   2.0268713 0.166666667
    ## OTU_56                   264    0.0000 0.00000000   0.9757790 0.047619048
    ## OTU_34                   781    0.0000 0.00000000   2.0268713 0.166666667
    ##         fisherAdjP     coeff      pvalues  adjPvalues  Kingdom
    ## OTU_10   1.0000000 -5.598472 2.453239e-03 0.008369873 Bacteria
    ## OTU_106  1.0000000 -3.848887 1.330121e-03 0.004821688 Bacteria
    ## OTU_18   1.0000000 -3.737111 2.428632e-04 0.001565119 Bacteria
    ## OTU_74   1.0000000 -3.606916 4.218682e-03 0.012234178 Bacteria
    ## OTU_28   1.0000000 -2.450339 5.757468e-04 0.002782776 Bacteria
    ## OTU_125  1.0000000 -2.381645 5.397998e-03 0.014908755 Bacteria
    ## OTU_71   1.0000000 -2.271126 5.448032e-04 0.002782776 Bacteria
    ## OTU_95   1.0000000 -2.081834 1.600566e-04 0.001565119 Bacteria
    ## OTU_11   1.0000000 -1.710174 1.284840e-03 0.004821688 Bacteria
    ## OTU_59   1.0000000 -1.692746 9.267702e-04 0.004134821 Bacteria
    ## OTU_70   1.0000000 -1.607536 9.350141e-03 0.023578617 Bacteria
    ## OTU_13   1.0000000 -1.606044 3.209666e-03 0.009797928 Bacteria
    ## OTU_42   1.0000000  1.861977 1.883088e-02 0.042007356 Bacteria
    ## OTU_73   1.0000000  3.104876 4.759308e-04 0.002760399 Bacteria
    ## OTU_27   1.0000000  3.178896 3.023464e-03 0.009742272 Bacteria
    ## OTU_82   1.0000000  3.580376 2.014315e-04 0.001565119 Bacteria
    ## OTU_66   1.0000000  3.789765 1.830101e-02 0.042007356 Bacteria
    ## OTU_38   1.0000000  4.005254 6.018323e-03 0.015866487 Bacteria
    ## OTU_9    1.0000000  4.134342 1.127444e-03 0.004670841 Bacteria
    ## OTU_24   1.0000000  4.981333 2.362143e-04 0.001565119 Bacteria
    ## OTU_54   1.0000000  5.135264 1.381184e-04 0.001565119 Bacteria
    ## OTU_68   0.4603175  5.547099 2.016248e-04 0.001565119 Bacteria
    ## OTU_64   1.0000000  6.090299 1.059880e-04 0.001565119 Bacteria
    ## OTU_56   0.5523810  6.575620 8.565233e-05 0.001565119 Bacteria
    ## OTU_34   1.0000000  6.644656 1.767143e-04 0.001565119 Bacteria
    ##                 Phylum               Class              Order
    ## OTU_10  Proteobacteria Gammaproteobacteria  Enterobacteriales
    ## OTU_106     Firmicutes          Clostridia      Clostridiales
    ## OTU_18      Firmicutes          Clostridia      Clostridiales
    ## OTU_74      Firmicutes     Erysipelotrichi Erysipelotrichales
    ## OTU_28      Firmicutes             Bacilli   Turicibacterales
    ## OTU_125     Firmicutes             Bacilli    Lactobacillales
    ## OTU_71      Firmicutes          Clostridia      Clostridiales
    ## OTU_95      Firmicutes          Clostridia      Clostridiales
    ## OTU_11      Firmicutes          Clostridia      Clostridiales
    ## OTU_59      Firmicutes          Clostridia      Clostridiales
    ## OTU_70      Firmicutes     Erysipelotrichi Erysipelotrichales
    ## OTU_13      Firmicutes          Clostridia      Clostridiales
    ## OTU_42   Bacteroidetes         Bacteroidia      Bacteroidales
    ## OTU_73      Firmicutes          Clostridia      Clostridiales
    ## OTU_27      Firmicutes          Clostridia      Clostridiales
    ## OTU_82      Firmicutes          Clostridia      Clostridiales
    ## OTU_66   Bacteroidetes         Bacteroidia      Bacteroidales
    ## OTU_38   Bacteroidetes         Bacteroidia      Bacteroidales
    ## OTU_9    Bacteroidetes         Bacteroidia      Bacteroidales
    ## OTU_24   Bacteroidetes         Bacteroidia      Bacteroidales
    ## OTU_54  Proteobacteria Gammaproteobacteria      Aeromonadales
    ## OTU_68  Actinobacteria      Coriobacteriia   Coriobacteriales
    ## OTU_64  Proteobacteria Gammaproteobacteria      Aeromonadales
    ## OTU_56   Bacteroidetes         Bacteroidia      Bacteroidales
    ## OTU_34   Bacteroidetes         Bacteroidia      Bacteroidales
    ##                        Family              Genus     Species
    ## OTU_10     Enterobacteriaceae               <NA>        <NA>
    ## OTU_106       Lachnospiraceae       Epulopiscium        <NA>
    ## OTU_18        Lachnospiraceae     [Ruminococcus]      gnavus
    ## OTU_74    Erysipelotrichaceae      [Eubacterium]    dolichum
    ## OTU_28      Turicibacteraceae       Turicibacter        <NA>
    ## OTU_125       Enterococcaceae       Enterococcus        <NA>
    ## OTU_71  Peptostreptococcaceae               <NA>        <NA>
    ## OTU_95        Lachnospiraceae              Dorea        <NA>
    ## OTU_11        Lachnospiraceae            Blautia    producta
    ## OTU_59        Lachnospiraceae     [Ruminococcus]        <NA>
    ## OTU_70    Erysipelotrichaceae               <NA>        <NA>
    ## OTU_13        Ruminococcaceae   Faecalibacterium prausnitzii
    ## OTU_42         Bacteroidaceae        Bacteroides coprophilus
    ## OTU_73        Ruminococcaceae       Oscillospira        <NA>
    ## OTU_27        Ruminococcaceae       Ruminococcus        <NA>
    ## OTU_82     [Mogibacteriaceae]               <NA>        <NA>
    ## OTU_66         Bacteroidaceae        Bacteroides   uniformis
    ## OTU_38     Porphyromonadaceae    Parabacteroides        <NA>
    ## OTU_9          Prevotellaceae         Prevotella       copri
    ## OTU_24                  S24-7               <NA>        <NA>
    ## OTU_54    Succinivibrionaceae               <NA>        <NA>
    ## OTU_68      Coriobacteriaceae      Adlercreutzia        <NA>
    ## OTU_64    Succinivibrionaceae Anaerobiospirillum        <NA>
    ## OTU_56                   <NA>               <NA>        <NA>
    ## OTU_34   [Paraprevotellaceae]               <NA>        <NA>

Now compare dog G vs. dog B (individual taxa)

    b = super.fitZig.kv(physeq = phy.temp,factor = "Dog",outDir = outDir,FileName =c("1_25FC_0.2_Dog_GvsB_OTUs"),
            heatmap.descriptor=c("tax_annot"), main=c("Dog G vs. B, OTUs"), subt=c("subt = FDR < 0.05,|coeff| >= 1.25, >20%+ in either group"), 
            ordered=TRUE, p=0.05, FC = 1.25, perc=0.2, extra.cols = c("Treatment"))

    ## [1] "0 of 10 samples were removed due to missing data"

    ## Default value being used.

    ## [1] "Dog will be modeled as a binary categorical predictor variable"

    ## Default value being used.

    ## it= 0, nll=18.24, log10(eps+1)=Inf, stillActive=138
    ## it= 1, nll=19.37, log10(eps+1)=Inf, stillActive=9
    ## it= 2, nll=19.41, log10(eps+1)=Inf, stillActive=3
    ## it= 3, nll=19.44, log10(eps+1)=Inf, stillActive=1
    ## it= 4, nll=19.44, log10(eps+1)=Inf, stillActive=1
    ## it= 5, nll=19.44, log10(eps+1)=Inf, stillActive=1
    ## it= 6, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it= 7, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it= 8, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it= 9, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=10, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=11, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=12, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=13, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=14, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=15, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=16, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=17, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=18, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=19, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=20, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=21, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=22, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=23, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=24, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=25, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=26, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=27, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=28, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=29, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=30, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=31, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=32, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=33, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=34, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=35, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=36, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=37, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=38, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=39, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=40, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=41, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=42, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=43, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=44, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=45, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=46, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=47, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=48, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=49, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=50, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=51, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=52, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=53, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=54, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=55, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=56, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=57, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=58, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=59, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=60, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=61, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=62, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=63, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=64, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=65, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=66, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=67, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=68, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=69, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=70, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=71, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=72, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=73, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=74, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=75, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=76, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=77, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=78, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=79, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=80, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=81, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=82, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=83, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=84, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=85, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=86, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=87, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=88, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=89, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=90, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=91, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=92, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=93, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=94, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=95, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=96, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=97, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=98, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## it=99, nll=19.43, log10(eps+1)=Inf, stillActive=1
    ## There were  44 OTUs significantly different between B vs. G that met 
    ##  threshold criteria of p 0.05 absolute FC 1.25 and percentage presence in at least one group of 20 % 
    ## [1] "writing results and model to file"
    ## [1] "/home/gerrit/scratch/amw/downstream/1_25FC_0.2_Dog_GvsB_OTUs_tax_annot.pdf"

![](README_files/figure-markdown_strict/unnamed-chunk-23-1.png)![](README_files/figure-markdown_strict/unnamed-chunk-23-2.png)

    ## [1] "making heatmap of results"

    b

    ##         percent_positive_group0 percent_positive_group1
    ## OTU_44                      100                     100
    ## OTU_45                      100                     100
    ## OTU_10                      100                     100
    ## OTU_71                      100                      80
    ## OTU_85                      100                      20
    ## OTU_57                      100                     100
    ## OTU_18                      100                     100
    ## OTU_70                      100                     100
    ## OTU_19                      100                     100
    ## OTU_20                      100                     100
    ## OTU_32                      100                     100
    ## OTU_28                      100                     100
    ## OTU_51                      100                     100
    ## OTU_11                      100                     100
    ## OTU_47                      100                     100
    ## OTU_13                      100                     100
    ## OTU_59                      100                     100
    ## OTU_92                       80                     100
    ## OTU_26                      100                     100
    ## OTU_104                     100                     100
    ## OTU_101                     100                     100
    ## OTU_141                     100                     100
    ## OTU_73                       80                     100
    ## OTU_52                       80                     100
    ## OTU_105                       0                     100
    ## OTU_174                      80                     100
    ## OTU_117                     100                     100
    ## OTU_9                       100                     100
    ## OTU_154                      40                      20
    ## OTU_82                       60                     100
    ## OTU_77                       60                     100
    ## OTU_24                      100                     100
    ## OTU_55                       20                     100
    ## OTU_96                       40                      40
    ## OTU_64                       40                     100
    ## OTU_68                        0                     100
    ## OTU_50                       40                     100
    ## OTU_66                       60                     100
    ## OTU_54                       60                     100
    ## OTU_78                       20                     100
    ## OTU_34                       40                     100
    ## OTU_56                       20                     100
    ## OTU_41                       40                     100
    ## OTU_116                     100                     100
    ##         +samples in group 0 +samples in group 1 mean_positive_group0
    ## OTU_44                    5                   5                  531
    ## OTU_45                    5                   5                  521
    ## OTU_10                    5                   5                 3645
    ## OTU_71                    5                   4                  136
    ## OTU_85                    5                   1                   36
    ## OTU_57                    5                   5                  152
    ## OTU_18                    5                   5                 1142
    ## OTU_70                    5                   5                   74
    ## OTU_19                    5                   5                 1320
    ## OTU_20                    5                   5                  583
    ## OTU_32                    5                   5                  103
    ## OTU_28                    5                   5                 1147
    ## OTU_51                    5                   5                  130
    ## OTU_11                    5                   5                 3377
    ## OTU_47                    5                   5                 1522
    ## OTU_13                    5                   5                 1822
    ## OTU_59                    5                   5                  105
    ## OTU_92                    4                   5                    4
    ## OTU_26                    5                   5                   98
    ## OTU_104                   5                   5                   19
    ## OTU_101                   5                   5                    6
    ## OTU_141                   5                   5                   32
    ## OTU_73                    4                   5                    6
    ## OTU_52                    4                   5                   32
    ## OTU_105                   0                   5                  NaN
    ## OTU_174                   4                   5                    5
    ## OTU_117                   5                   5                   14
    ## OTU_9                     5                   5                 1510
    ## OTU_154                   2                   1                    2
    ## OTU_82                    3                   5                    7
    ## OTU_77                    3                   5                    3
    ## OTU_24                    5                   5                   20
    ## OTU_55                    1                   5                    3
    ## OTU_96                    2                   2                   22
    ## OTU_64                    2                   5                    2
    ## OTU_68                    0                   5                  NaN
    ## OTU_50                    2                   5                    2
    ## OTU_66                    3                   5                    4
    ## OTU_54                    3                   5                    2
    ## OTU_78                    1                   5                    1
    ## OTU_34                    2                   5                   10
    ## OTU_56                    1                   5                    4
    ## OTU_41                    2                   5                    4
    ## OTU_116                   5                   5                   23
    ##         mean_positive_group1 oddsRatio      lower       upper     fisherP
    ## OTU_44                     9  0.000000 0.00000000         Inf 1.000000000
    ## OTU_45                    13  0.000000 0.00000000         Inf 1.000000000
    ## OTU_10                    70  0.000000 0.00000000         Inf 1.000000000
    ## OTU_71                     4       Inf 0.02564066         Inf 1.000000000
    ## OTU_85                     2       Inf 1.02482226         Inf 0.047619048
    ## OTU_57                    17  0.000000 0.00000000         Inf 1.000000000
    ## OTU_18                    85  0.000000 0.00000000         Inf 1.000000000
    ## OTU_70                     8  0.000000 0.00000000         Inf 1.000000000
    ## OTU_19                   115  0.000000 0.00000000         Inf 1.000000000
    ## OTU_20                   245  0.000000 0.00000000         Inf 1.000000000
    ## OTU_32                    77  0.000000 0.00000000         Inf 1.000000000
    ## OTU_28                   241  0.000000 0.00000000         Inf 1.000000000
    ## OTU_51                    50  0.000000 0.00000000         Inf 1.000000000
    ## OTU_11                   806  0.000000 0.00000000         Inf 1.000000000
    ## OTU_47                   554  0.000000 0.00000000         Inf 1.000000000
    ## OTU_13                   588  0.000000 0.00000000         Inf 1.000000000
    ## OTU_59                    32  0.000000 0.00000000         Inf 1.000000000
    ## OTU_92                    26  0.000000 0.00000000  39.0005500 1.000000000
    ## OTU_26                  1026  0.000000 0.00000000         Inf 1.000000000
    ## OTU_104                   86  0.000000 0.00000000         Inf 1.000000000
    ## OTU_101                   32  0.000000 0.00000000         Inf 1.000000000
    ## OTU_141                  147  0.000000 0.00000000         Inf 1.000000000
    ## OTU_73                    48  0.000000 0.00000000  39.0005500 1.000000000
    ## OTU_52                   194  0.000000 0.00000000  39.0005500 1.000000000
    ## OTU_105                    8  0.000000 0.00000000   0.4353226 0.007936508
    ## OTU_174                   20  0.000000 0.00000000  39.0005500 1.000000000
    ## OTU_117                  152  0.000000 0.00000000         Inf 1.000000000
    ## OTU_9                  15388  0.000000 0.00000000         Inf 1.000000000
    ## OTU_154                   14  2.414224 0.08474680 195.6529809 1.000000000
    ## OTU_82                    72  0.000000 0.00000000   5.1183766 0.444444444
    ## OTU_77                    45  0.000000 0.00000000   5.1183766 0.444444444
    ## OTU_24                   455  0.000000 0.00000000         Inf 1.000000000
    ## OTU_55                    93  0.000000 0.00000000   0.9757790 0.047619048
    ## OTU_96                    19  1.000000 0.04224561  23.6710987 1.000000000
    ## OTU_64                   138  0.000000 0.00000000   2.0268713 0.166666667
    ## OTU_68                    52  0.000000 0.00000000   0.4353226 0.007936508
    ## OTU_50                    88  0.000000 0.00000000   2.0268713 0.166666667
    ## OTU_66                    88  0.000000 0.00000000   5.1183766 0.444444444
    ## OTU_54                   104  0.000000 0.00000000   5.1183766 0.444444444
    ## OTU_78                    52  0.000000 0.00000000   0.9757790 0.047619048
    ## OTU_34                   781  0.000000 0.00000000   2.0268713 0.166666667
    ## OTU_56                   264  0.000000 0.00000000   0.9757790 0.047619048
    ## OTU_41                   635  0.000000 0.00000000   2.0268713 0.166666667
    ## OTU_116                 6243  0.000000 0.00000000         Inf 1.000000000
    ##         fisherAdjP     coeff      pvalues   adjPvalues  Kingdom
    ## OTU_44   1.0000000 -6.452149 1.278213e-04 2.188940e-03 Bacteria
    ## OTU_45   1.0000000 -6.258280 5.161248e-08 7.070910e-06 Bacteria
    ## OTU_10   1.0000000 -4.822274 6.342219e-03 2.473993e-02 Bacteria
    ## OTU_71   1.0000000 -4.161072 6.863480e-04 5.531158e-03 Bacteria
    ## OTU_85   0.4693878 -3.675842 7.287746e-04 5.546784e-03 Bacteria
    ## OTU_57   1.0000000 -3.166523 1.690626e-03 1.007025e-02 Bacteria
    ## OTU_18   1.0000000 -3.153951 1.122893e-03 7.325541e-03 Bacteria
    ## OTU_70   1.0000000 -3.140823 2.647313e-03 1.450727e-02 Bacteria
    ## OTU_19   1.0000000 -3.039097 9.019078e-05 2.028373e-03 Bacteria
    ## OTU_20   1.0000000 -2.731507 3.780328e-03 1.803105e-02 Bacteria
    ## OTU_32   1.0000000 -2.704828 8.567859e-03 3.088939e-02 Bacteria
    ## OTU_28   1.0000000 -2.577207 1.431669e-03 8.915395e-03 Bacteria
    ## OTU_51   1.0000000 -2.065738 9.319273e-04 6.719686e-03 Bacteria
    ## OTU_11   1.0000000 -1.982204 9.617979e-03 3.294158e-02 Bacteria
    ## OTU_47   1.0000000 -1.949285 1.853111e-03 1.057817e-02 Bacteria
    ## OTU_13   1.0000000 -1.801759 6.011705e-03 2.422363e-02 Bacteria
    ## OTU_59   1.0000000 -1.698921 5.628693e-03 2.336760e-02 Bacteria
    ## OTU_92   1.0000000  1.649828 9.601959e-03 3.294158e-02 Bacteria
    ## OTU_26   1.0000000  1.943959 1.210695e-02 4.045492e-02 Bacteria
    ## OTU_104  1.0000000  2.061113 1.388219e-02 4.322410e-02 Bacteria
    ## OTU_101  1.0000000  2.212144 5.331985e-03 2.282756e-02 Bacteria
    ## OTU_141  1.0000000  2.297183 6.501002e-03 2.473993e-02 Bacteria
    ## OTU_73   1.0000000  2.469124 3.725810e-03 1.803105e-02 Bacteria
    ## OTU_52   1.0000000  2.609183 1.339480e-02 4.267645e-02 Bacteria
    ## OTU_105  0.2738095  2.807601 6.746812e-03 2.498144e-02 Bacteria
    ## OTU_174  1.0000000  2.849916 3.816792e-03 1.803105e-02 Bacteria
    ## OTU_117  1.0000000  3.073081 1.006966e-03 6.897716e-03 Bacteria
    ## OTU_9    1.0000000  3.318921 4.471045e-03 2.041777e-02 Bacteria
    ## OTU_154  1.0000000  3.334589 1.318015e-02 4.267645e-02 Bacteria
    ## OTU_82   1.0000000  3.479223 3.618063e-04 3.831248e-03 Bacteria
    ## OTU_77   1.0000000  3.502925 5.132801e-04 4.687958e-03 Bacteria
    ## OTU_24   1.0000000  4.110471 6.253611e-04 5.354654e-03 Bacteria
    ## OTU_55   0.4693878  4.564643 3.424357e-04 3.831248e-03 Bacteria
    ## OTU_96   1.0000000  4.728294 4.799093e-03 2.120889e-02 Bacteria
    ## OTU_64   1.0000000  4.929019 3.385932e-04 3.831248e-03 Bacteria
    ## OTU_68   0.2738095  4.966152 3.635491e-04 3.831248e-03 Bacteria
    ## OTU_50   1.0000000  5.017722 2.461335e-04 3.746699e-03 Bacteria
    ## OTU_66   1.0000000  5.047824 3.648418e-03 1.803105e-02 Bacteria
    ## OTU_54   1.0000000  5.277387 5.616288e-05 1.538863e-03 Bacteria
    ## OTU_78   0.4693878  5.486823 4.784650e-05 1.538863e-03 Bacteria
    ## OTU_34   1.0000000  5.766382 1.036395e-04 2.028373e-03 Bacteria
    ## OTU_56   0.4693878  5.841937 4.867022e-05 1.538863e-03 Bacteria
    ## OTU_41   1.0000000  6.100335 4.439625e-04 4.344490e-03 Bacteria
    ## OTU_116  1.0000000  8.970073 6.057725e-07 4.149542e-05 Bacteria
    ##                 Phylum               Class              Order
    ## OTU_44      Firmicutes          Clostridia      Clostridiales
    ## OTU_45      Firmicutes          Clostridia      Clostridiales
    ## OTU_10  Proteobacteria Gammaproteobacteria  Enterobacteriales
    ## OTU_71      Firmicutes          Clostridia      Clostridiales
    ## OTU_85      Firmicutes     Erysipelotrichi Erysipelotrichales
    ## OTU_57      Firmicutes          Clostridia      Clostridiales
    ## OTU_18      Firmicutes          Clostridia      Clostridiales
    ## OTU_70      Firmicutes     Erysipelotrichi Erysipelotrichales
    ## OTU_19      Firmicutes          Clostridia      Clostridiales
    ## OTU_20   Bacteroidetes         Bacteroidia      Bacteroidales
    ## OTU_32   Bacteroidetes         Bacteroidia      Bacteroidales
    ## OTU_28      Firmicutes             Bacilli   Turicibacterales
    ## OTU_51      Firmicutes          Clostridia      Clostridiales
    ## OTU_11      Firmicutes          Clostridia      Clostridiales
    ## OTU_47      Firmicutes          Clostridia      Clostridiales
    ## OTU_13      Firmicutes          Clostridia      Clostridiales
    ## OTU_59      Firmicutes          Clostridia      Clostridiales
    ## OTU_92      Firmicutes          Clostridia      Clostridiales
    ## OTU_26      Firmicutes          Clostridia      Clostridiales
    ## OTU_104   Fusobacteria       Fusobacteriia    Fusobacteriales
    ## OTU_101     Firmicutes          Clostridia      Clostridiales
    ## OTU_141   Fusobacteria       Fusobacteriia    Fusobacteriales
    ## OTU_73      Firmicutes          Clostridia      Clostridiales
    ## OTU_52   Bacteroidetes         Bacteroidia      Bacteroidales
    ## OTU_105     Firmicutes     Erysipelotrichi Erysipelotrichales
    ## OTU_174     Firmicutes          Clostridia      Clostridiales
    ## OTU_117  Bacteroidetes         Bacteroidia      Bacteroidales
    ## OTU_9    Bacteroidetes         Bacteroidia      Bacteroidales
    ## OTU_154     Firmicutes          Clostridia      Clostridiales
    ## OTU_82      Firmicutes          Clostridia      Clostridiales
    ## OTU_77      Firmicutes          Clostridia      Clostridiales
    ## OTU_24   Bacteroidetes         Bacteroidia      Bacteroidales
    ## OTU_55  Proteobacteria Gammaproteobacteria      Aeromonadales
    ## OTU_96      Firmicutes          Clostridia      Clostridiales
    ## OTU_64  Proteobacteria Gammaproteobacteria      Aeromonadales
    ## OTU_68  Actinobacteria      Coriobacteriia   Coriobacteriales
    ## OTU_50   Bacteroidetes         Bacteroidia      Bacteroidales
    ## OTU_66   Bacteroidetes         Bacteroidia      Bacteroidales
    ## OTU_54  Proteobacteria Gammaproteobacteria      Aeromonadales
    ## OTU_78      Firmicutes          Clostridia      Clostridiales
    ## OTU_34   Bacteroidetes         Bacteroidia      Bacteroidales
    ## OTU_56   Bacteroidetes         Bacteroidia      Bacteroidales
    ## OTU_41   Bacteroidetes         Bacteroidia      Bacteroidales
    ## OTU_116  Bacteroidetes         Bacteroidia      Bacteroidales
    ##                        Family              Genus     Species
    ## OTU_44        Lachnospiraceae              Dorea            
    ## OTU_45        Lachnospiraceae              Dorea            
    ## OTU_10     Enterobacteriaceae                               
    ## OTU_71  Peptostreptococcaceae                               
    ## OTU_85    Erysipelotrichaceae                               
    ## OTU_57         Clostridiaceae                               
    ## OTU_18        Lachnospiraceae     [Ruminococcus]      gnavus
    ## OTU_70    Erysipelotrichaceae                               
    ## OTU_19        Lachnospiraceae                               
    ## OTU_20         Bacteroidaceae        Bacteroides            
    ## OTU_32         Bacteroidaceae        Bacteroides            
    ## OTU_28      Turicibacteraceae       Turicibacter            
    ## OTU_51        Lachnospiraceae                               
    ## OTU_11        Lachnospiraceae            Blautia    producta
    ## OTU_47        Lachnospiraceae               <NA>        <NA>
    ## OTU_13        Ruminococcaceae   Faecalibacterium prausnitzii
    ## OTU_59        Lachnospiraceae     [Ruminococcus]            
    ## OTU_92        Ruminococcaceae       Oscillospira            
    ## OTU_26        Lachnospiraceae                               
    ## OTU_104      Fusobacteriaceae      Fusobacterium            
    ## OTU_101       Ruminococcaceae       Oscillospira            
    ## OTU_141      Fusobacteriaceae      Fusobacterium            
    ## OTU_73        Ruminococcaceae       Oscillospira            
    ## OTU_52         Bacteroidaceae        Bacteroides            
    ## OTU_105   Erysipelotrichaceae                               
    ## OTU_174       Lachnospiraceae              Dorea            
    ## OTU_117        Prevotellaceae         Prevotella       copri
    ## OTU_9          Prevotellaceae         Prevotella       copri
    ## OTU_154       Ruminococcaceae                               
    ## OTU_82     [Mogibacteriaceae]                               
    ## OTU_77        Lachnospiraceae                               
    ## OTU_24                  S24-7                               
    ## OTU_55    Succinivibrionaceae Anaerobiospirillum            
    ## OTU_96        Veillonellaceae        Megasphaera            
    ## OTU_64    Succinivibrionaceae Anaerobiospirillum            
    ## OTU_68      Coriobacteriaceae      Adlercreutzia            
    ## OTU_50     Porphyromonadaceae    Parabacteroides            
    ## OTU_66         Bacteroidaceae        Bacteroides   uniformis
    ## OTU_54    Succinivibrionaceae                               
    ## OTU_78                                                      
    ## OTU_34   [Paraprevotellaceae]               <NA>        <NA>
    ## OTU_56                                                      
    ## OTU_41   [Paraprevotellaceae]       [Prevotella]            
    ## OTU_116  [Paraprevotellaceae]       [Prevotella]

    sessionInfo()

    ## R version 3.3.2 (2016-10-31)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 14.04.5 LTS
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_ZA.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_ZA.UTF-8        LC_COLLATE=en_ZA.UTF-8    
    ##  [5] LC_MONETARY=en_ZA.UTF-8    LC_MESSAGES=en_ZA.UTF-8   
    ##  [7] LC_PAPER=en_ZA.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_ZA.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] ROCR_1.0-7           gplots_3.0.1         dplyr_0.5.0         
    ##  [4] randomForest_4.6-12  metagenomeSeq_1.16.0 RColorBrewer_1.1-2  
    ##  [7] glmnet_2.0-5         foreach_1.4.3        Matrix_1.2-8        
    ## [10] limma_3.30.13        fifer_1.1            MASS_7.3-45         
    ## [13] matrixStats_0.51.0   psych_1.7.3.21       corrplot_0.77       
    ## [16] vegan_2.4-2          lattice_0.20-35      permute_0.9-4       
    ## [19] NMF_0.20.6           Biobase_2.34.0       BiocGenerics_0.20.0 
    ## [22] cluster_2.0.6        rngtools_1.2.4       pkgmaker_0.22       
    ## [25] registry_0.3         dunn.test_1.3.3      gridExtra_2.2.1     
    ## [28] ggplot2_2.2.1        phyloseq_1.19.1     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] nlme_3.1-131          bitops_1.0-6          doParallel_1.0.10    
    ##  [4] rprojroot_1.2         tools_3.3.2           backports_1.0.5      
    ##  [7] R6_2.2.0              KernSmooth_2.23-15    rpart_4.1-10         
    ## [10] DBI_0.6               Hmisc_4.0-2           lazyeval_0.2.0       
    ## [13] mgcv_1.8-17           colorspace_1.3-2      ade4_1.7-6           
    ## [16] nnet_7.3-12           mnormt_1.5-5          compiler_3.3.2       
    ## [19] htmlTable_1.9         sandwich_2.3-4        labeling_0.3         
    ## [22] caTools_1.17.1        scales_0.4.1          checkmate_1.8.2      
    ## [25] mvtnorm_1.0-6         stringr_1.2.0         digest_0.6.12        
    ## [28] foreign_0.8-67        rmarkdown_1.4         XVector_0.14.1       
    ## [31] base64enc_0.1-3       htmltools_0.3.5       plotrix_3.6-4        
    ## [34] maps_3.1.1            htmlwidgets_0.8       zoo_1.7-14           
    ## [37] jsonlite_1.3          gtools_3.5.0          acepack_1.4.1        
    ## [40] magrittr_1.5          modeltools_0.2-21     Formula_1.2-1        
    ## [43] biomformat_1.2.0      Rcpp_0.12.10          munsell_0.4.3        
    ## [46] S4Vectors_0.12.2      ape_4.1               stringi_1.1.3        
    ## [49] multcomp_1.4-6        yaml_2.1.14           zlibbioc_1.20.0      
    ## [52] rhdf5_2.18.0          plyr_1.8.4            grid_3.3.2           
    ## [55] strucchange_1.5-1     gdata_2.17.0          Biostrings_2.42.1    
    ## [58] splines_3.3.2         multtest_2.30.0       knitr_1.15.1         
    ## [61] igraph_1.0.1          party_1.2-2           reshape2_1.4.2       
    ## [64] codetools_0.2-15      stats4_3.3.2          evaluate_0.10        
    ## [67] latticeExtra_0.6-28   data.table_1.10.4     spam_1.4-0           
    ## [70] gtable_0.2.0          assertthat_0.1        gridBase_0.4-7       
    ## [73] coin_1.1-3            xtable_1.8-2          survival_2.41-2      
    ## [76] randomForestSRC_2.4.2 tibble_1.2            iterators_1.0.8      
    ## [79] IRanges_2.8.2         fields_8.10           TH.data_1.0-8

Packages req required run this tutorial
---------------------------------------

-   phyloseq
-   ggplot2
-   gridExtra
-   dunn.test
-   NMF
-   vegan
-   corrplot
-   psych
-   matrixStats
-   fifer
-   metagenomeSeq
-   randomForest
-   dplyr
-   ROCR
