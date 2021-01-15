library(phyloseq)
library(NMF)
library(vegan)
library(ggplot2)
library(corrplot)#for cor.mtest
library(psych)#corr.test
library(matrixStats)#rowSds
library(metagenomeSeq)#differential abundance testing
library(randomForest)
library(dplyr)
library(ROCR)
library(caret)
library(gridExtra)#for grid.arrange()
nmf.options(grid.patch=TRUE)#set to avoid blank first pdf page being created

#DEFINE CUSTOM COLOR PALETTE WHERE COLORS ARE EASIER TO DISTINGUISH FROM ONE ANOTHER
myPalette <- c('#89C5DA', "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861")

#???????????????????????????????????????????????????????????
#SUMMARY OF FUNCTIONS AVAILABLE:
#1. make_metagenomeSeq() - converts phyloseq to metagenomeSeq object
#2. tax.lab() - GET LOWEST LEVEL OF TAXONOMIC ANNOTATION AVAILABLE TO SUBSTITUTE OTUS FOR PLOTTING PURPOSES.
#3. heatmap.k() - Generic heatmap function for phyloseq objects using package NMF
#4. bar.plots() - pretty generic barplot function built on phyloseq's plot_bar()
#5. tax_glom.kv() - GET LOWEST LEVEL OF TAXONOMIC ANNOTATION AVAILABLE AND MERGE COUNTS AT THIS LEVEL 
#6. super.fitZig.kv() - FUNCTION BUILT AROUND METAGENOMESEQ'S FITZIG() AND MRFULLTABLE() FUNCTIONS INCLUDING CUSTOM FILTERING AND HEATMAP OF SIGNIFICANT RESULTS (currently only setup for 2-class comparisons)
#7. corr.k(): FUNCTION TO PERFORM PAIRWISE CORRELATIONS, TEST FOR SIGNIFICANCE AND PLOT CORRELLOGRAM
#8. MC.summary(): take the output from PICRUSt's metagenome_contributions.py, together with taxonomic annotation for the OTUs included in 
#9. #RF.k(): performs random forests analysis on the otu table of a supplied phyloseq object.
#NB: this function is only setup for two-class categorical response variables NOT regression on continuous response variables

#Note: if getting error with NMF aheatmap() see https://github.com/renozao/NMF/issues/96 for troubleshooting
#???????????????????????????????????????????????????????????

#**********************************
#make_metagenomeSeq: THIS METAGENOME SEQ FUNCTION WAS COPIED AND MODIFIED FROM CODE USED FOR THE MCMURDIE 2014 PAPER - http://joey711.github.io/waste-not-supplemental/simulation-differential-abundance/simulation-differential-abundance-server.html
#**********************************
# Function to convert from phyloseq object to metagenomeSeq object
make_metagenomeSeq = function(physeq) {
  require("metagenomeSeq")
  require("phyloseq")
  # Enforce orientation
  if (!taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }
  OTU = as(otu_table(physeq), "matrix")
  # Convert sample_data to AnnotatedDataFrame
  ADF = AnnotatedDataFrame(data.frame(sample_data(physeq)))
  # define dummy 'feature' data for OTUs, using their name Helps with
  # extraction and relating to taxonomy later on.
  TDF = AnnotatedDataFrame(data.frame(OTUname = taxa_names(physeq), row.names = taxa_names(physeq)))
  # Create the metagenomeSeq object
  MGS = newMRexperiment(counts = OTU, phenoData = ADF, featureData = TDF)
  # Trigger metagenomeSeq to calculate its Cumulative Sum scaling factor.
  MGS = cumNorm(MGS)
  return(MGS)
}

#**********************************
#tax.lab() - GET LOWEST LEVEL OF TAXONOMIC ANNOTATION AVAILABLE TO SUBSTITUTE OTUS FOR PLOTTING PURPOSES.
#**********************************
#ARGUMENTS:
#otus = otu table of interest (often subset)
#physeq = phyloseq object of interst (for tax annotation)
#labrow: TRUE=annotate rows with taxonomy info, FALSE=annotate rows with OTU info.
#merged: if OTUs have been merged, level of merge e.g. "Genus", else FALSE
tax.lab <- function(otus, physeq, labrow=TRUE,merged=FALSE){
  #First delete columns with only NAs (generated for merged tax tables e.g. at Family level will have NA Genus, Species cols)
  temp_tax <- tax_table(physeq)
  temp_tax <- temp_tax[, colSums(is.na(temp_tax)) != nrow(temp_tax)]
  tax <- data.frame(temp_tax)
  tax <- tax[c(rownames(otus)),]
  tax[is.na(tax)] <- ""
  if(merged==FALSE){
    tax$unique <- rep("",dim(tax)[1])#make unique names
    for(i in 1:dim(tax)[1]){
      x = 7#column 7 = Species
      tax$unique[i] <- ifelse(tax[i,x]==""|tax[i,x]=="NA"|is.na(tax[i,x])==TRUE,
                              ifelse(tax[i,x-1]==""|tax[i,x-1]=="NA"|is.na(tax[i,x-1])==TRUE,
                                     ifelse(tax[i,x-2]==""|tax[i,x-2]=="NA"|is.na(tax[i,x-2])==TRUE,
                                            ifelse(tax[i,x-3]==""|tax[i,x-3]=="NA"|is.na(tax[i,x-3])==TRUE,#phylum
                                                   ifelse(tax[i,x-4]==""|tax[i,x-4]=="NA"|is.na(tax[i,x-4])==TRUE, as.character(tax[i,x-5]),
                                                          as.character(tax[i,x-4])),#class
                                                   as.character(tax[i,x-3])),#order
                                            as.character(tax[i,x-2])),#family
                                     as.character(tax[i,x-1])),#genus
                              as.character(tax[i,x]))#species
    }
    #Custom cleanup - species:
    #For OTUs with species-level annotation: get the first character of the Genus and append to species.
    #first remove all "[]"brackets
    tax$Genus = sub("\\[","",tax$Genus)
    tax$Genus = sub("\\]","",tax$Genus)
    for(i in 1:dim(tax)[1]){
      if(is.na(tax$Species[i])==FALSE & tax$Species[i]!="" & tax$Species[i]!="NA"){#get OTUs with species-level annotation
        get = substr(tax$Genus[i],1,1)
        tax$unique[i] = sub(tax$unique[i],paste0(get,".",tax$unique[i]),tax$unique[i])
      }	
    }
    tax$unique <- make.unique(tax$unique)
    
    if(labrow==TRUE){
      labs = tax$unique
    }else if(labrow==FALSE){
      labs = rownames(otus)
    }
  }else if(length(merged)>0){
    if(labrow==TRUE){
      labs=tax[,merged]
    }else if(labrow==FALSE){
      labs = rownames(otus)
    }
    
  }
  return(labs)
}

#**********************************
#heatmap.k() Generic heatmap function for phyloseq 
#objects using package NMF
#**********************************
#ARGUMENTS:
#physeq = phyloseq object of interst (for tax annotation) - only applicable if you're plotting OTUs.
#-
#plot.otus = are you plotting otus, true or false, default = TRUE
#-
#data.subset - either a dataframe with rownames = the otus you want to plot (typically the significant results table from metagenomeseq)
#OR a dataframe with some other data you want to plot e.g. cytokines, flow cyt, gene expression
#-
#subset.samples = names of samples to include in heatmap, default = NULL, i.e. use all samples
#annot.cols = columns of interest in sample_data(physeq) for annotation track
#colours = typical NMF color specification for columns - see NMF reference manual
#order.by = column number by which to order columns (variable of interest), default = NULL (i.e. perform hierarchical clustering)
#main = Title of plot (analysis details)
#subt = subtitle of plot (threshold info)
#filename = filename
#Colv = should columns be clustered, default = no
#labrow: TRUE=annotate rows with taxonomy info, FALSE=annotate rows with OTU info. - only applicable if you're plotting OTUs
#if OTUs have been merged, TRUE, else FALSE - only applicable if you're plotting OTUs
#merged = if the data has been merged at a given taxonomic level - supply column name e.g. "Genus"
#distfun = distance function used, default 'euclidian'
#hclustfun = clustering function used, default 'average'
#scale: see aheatmap() from NMF package for details. Should heatmap rows/columns be scaled?

heatmap.k = function(physeq=NULL, plot.otus = TRUE,data.subset=NULL,subset.samples = NULL, 
                     annot.cols = NULL,colours = NULL, order.by = NULL,main=NULL, subt=NULL,filename,
                     Colv = NULL,rows.sortby=NULL,labrow = FALSE, merged = FALSE, distfun = 'euclidean',hclustfun = 'average', cexCol=1, cexRow=1, scale = "none"){
  if(plot.otus==TRUE){
    #If plotting OTUs - get OTUs
    otus <- otu_table(physeq)
    if(length(subset.samples)==0){
      print("including all samples")#do nothing
    }else if(length(subset.samples) >0){#if there was a subset specified:
      physeq <- prune_samples(sample_names(physeq)%in% subset.samples,physeq)
      otus <- otus[,colnames(otus) %in% subset.samples]
    }
    if(length(data.subset)==0){
      print("including all otus")#do nothing
    }else if(length(data.subset) >0){#if an otus subset was specified
      otus <- otus[rownames(otus) %in% rownames(data.subset),]
    }
    #log2 transform otus for better heatmap visualisation colour scale
    zs <- which(otus==0)#find zero entries 
    if(length(zs)>0){#if there are zero values, transform
      otus[zs] <- 1 #set zs to 1 pseudocount
      otus <- log2(otus)
    } else{otus <- log2(otus)}
    plot <- otus#data to be plotted
    #--
  }else if(plot.otus==FALSE){
    plot = data.subset#this is for non-otu data
    if(length(colnames(plot) %in% sample_names(physeq))>0){
      #make sure the table is in the correct orientation
      
    }else{plot = t(plot)}#transform so that sample names are column names not rownames
    #make sure to subset the data appropriately to include only samples in data.subset
    sample_data(physeq) <- sample_data(physeq)[rownames(sample_data(physeq)) %in% colnames(plot),]
    
  }
  #get annotation
  pheno <- sample_data(physeq)
  a.track <- pheno[,annot.cols]#select columns of interst
  if(length(order.by)==0){
    #do nothing
  }else {a.track <- a.track[order(a.track[,order.by][[1]], decreasing = TRUE),]}#sort by variable of interest
  #make sure the samples used in the data to plot match that of the annotation track 
  plot = plot[,rownames(a.track)]
  #sort ROWS by column of interest? (e.g. if you want to sort by p-value or coefficient size)
  if(length(rows.sortby)!=0){
    #this is only for results from super.fitZig.kv and makes use of the column specified in 'rows.sortby'
    if(length(which(colnames(data.subset)==rows.sortby))==1){#
      ids.order <- rownames(data.subset[order(abs(data.subset[,rows.sortby])),])
      plot <- plot[ids.order,]
      rows.sortby = NA	
    }
  }
  #-------------------
  #row tax annotation: if plotting otus
  if(plot.otus==TRUE){
    labs <- tax.lab(physeq=physeq,otus=plot,labrow, merged = merged)
  }else{labs = rownames(plot)}
  #Write to file
  nmf.options(grid.patch=TRUE)#set to avoid blank first pdf page being created
  pdf(filename)
  a <- aheatmap(plot, main = main,			
                sub = subt,
                scale = scale,
                Colv = Colv,
                Rowv = rows.sortby,
                annCol = a.track,
                annColors = colours,
                labRow = labs,
                info = "distance",
                distfun = distfun,
                hclustfun=hclustfun,
                fontsize=10,
                cexCol=cexCol,
                cexRow=cexRow)
  dev.off()
  #redo a, otherwise not returned
  a <- aheatmap(plot, main = main,
                sub = subt,
                Colv = Colv,
                scale = scale,
                Rowv = rows.sortby,
                annCol = a.track,
                annColors = colours,
                labRow = labs,
                info = "distance",
                distfun = distfun,
                hclustfun=hclustfun,
                cexCol=cexCol,
                cexRow=cexRow)
  return(a)	
}

#**********************************
#bar.plots() - generic barplot function build on phyloseq plot_bar(): 
#NB number of taxa that can be displayed currently limited to 26 (number defined in myPalette at start of script)
#**********************************
#This function was modified from the phyloseq plot_bar() function where ggplot2's geom_bar no longers sorted stacked bars by abundance.
#Don't call plot_bar_fix directly, it gets called by bar.plots()
#---
#ARGUMENTS for bar.plots()
#physeq: phyloseq object, standardized counts
#cat: category of interest
#level: level at which to perform taxa merging
#perc: minimum fraction of samples which should be positive for a given taxa (stringent filtering (perc=0.5) used to minimize number of taxa in figure legend, and only display most abundant taxa)
#subset.otus = any dataframe for which the rownames will be used to subset otu.table (stats.sig)
#subset.samples = names of samples to include in barplot
#count = minimum per taxon count summed across all samples
#order.bars = option to order barplots  - sometimes relavant if more than 2 categories (e.g. order= c("classA","classB","classC")
#x.axis =x axis labels, 
#y.axis = y axis labels, 
#filen = any other filename additions such as project name or date

plot_bar_fix <- function (physeq, x = "Sample", y = "Abundance", fill = level, 
                          title = NULL, facet_grid = NULL) 
{
  mdf = psmelt(physeq)
  p=ggplot(mdf[order(-mdf[,y]),],aes_string(x= x, y=y, fill = fill, group = y)) + geom_bar(stat='identity')
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

bar.plots <- function(physeq,subset.otus  = NULL, subset.samples = NULL, count,cat, 
                      level, perc, order.bars = NULL,x.axis=NULL, y.axis=NULL, merge.samples = FALSE,filen = "",outDir){
  if(length(subset.samples)==0){
    
  }else if (length(subset.samples) >0){#if there was a subset specified:
    physeq <- prune_samples(sample_names(physeq)%in% subset.samples,physeq)
  }
  if(length(subset.otus)==0){
    
  }else if(dim(subset.otus)[1] >0){#if an otus subset was specified
    physeq <- prune_taxa(rownames(tax_table(physeq)) %in% rownames(subset.otus),physeq)
  }
  #MERGE OTUS AT USER-SPECIFIED LEVEL
  Ns=table(sample_data(physeq)[,cat])#get N's for each category to be listed in plot title
  merged <- tax_glom(physeq, level)
  #filter taxa using the same filter used for diff. abundance testing
  merged = filter_taxa(merged, function(x) sum(x > count) >= (perc*length(x)), TRUE)#filter taxa for barplot display
  if(merge.samples==FALSE){
    new <- transform_sample_counts(merged, function(x) 100 * x/sum(x))#convert to % abundance
    
  }else{
    new <- merge_samples(merged,cat, fun = "median")#Merge by category
    new <- transform_sample_counts(new, function(x) 100 * x/sum(x))#convert to % abundance
    sample_data(new)[,cat] <- rownames(sample_data(new)) #merge_samples also "merges the sample table, replace relevant cat column with rownames
  }
 
  #make vector with number of samples per group to include in title
  if(length(order.bars)==0){#If and order is specified we also need to change the order of the group counts displayed in for.title
    for.title = c()
    for(i in 1: length(Ns)){
      for.title <- c(Ns[[i]],for.title)
    }
  }else{
    for.title = c()
    Ns = Ns[order.bars]
    for(i in 1: length(Ns)){
      for.title <- c(Ns[[i]],for.title)
    }
  }	
  for.title <- paste(for.title, collapse = ",")#now paste as one string for use in title
  p2 <- plot_bar_fix(physeq=new,x = cat, fill=level ,title = paste0(level,"-level 16S ",x.axis,". Ns = ",for.title))+ coord_flip() + 
    ylab("Percentage of Sequences")+theme(axis.text=element_text(size=14),
                                          axis.title=element_text(size=14),legend.title=element_text(size=16),
                                          legend.text=element_text(size=14),
                                          plot.title=element_text(size=14, face="bold"))+scale_x_discrete(limits=c(order.bars))+scale_fill_manual(values=myPalette)+theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  p2 = p2+labs(x.axis=x.axis, y.axis=y.axis)
  pdf(paste0(outDir,"/",filen,"",level,"_abundance_by_",x.axis,perc,"_",count,".pdf"))
  par(cex.lab = 2)
  par(cex.axis = 2)
  grid.arrange(p2, ncol=1)	
  dev.off()
  return(p2)
  
}

#**********************************
##tax_glom.kv - combination of phyloseq's tax_glom() and my tax.lab() function
#GET LOWEST LEVEL OF TAXONOMIC ANNOTATION AVAILABLE AND MERGE COUNTS AT THIS LEVEL (i.e. I have 3 otus of L.iners and another 3 of Lactobacillus (unknown species) -->
#e.g. merge L.iners at species level and Lactobacillus (no species annot) at genus level.
#modified from phyloseq function tax_glom
#NB: this works for when we have the standard 7 ranks
#colnames(tax_table(physeq))
#[1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"
#**********************************
#ARGUMENTS:
#otus = otu table of interest, if not specified, defaults to OTU table of the phyloseq obj. specified.
#physeq = phyloseq object with taxonomic annotation
tax_glom.kv <- function (physeq) 
{
  otus = otu_table(physeq)
  if (is.null(access(physeq, "tax_table"))) {
    stop("The tax_glom.kv() function requires that physeq contain a taxonomyTable")
  }
  #check if there is a sample_data slot - function merge_phyloseq() downstream requires a sample_data slot (bug?)
  remove.SD = FALSE
  remove.SD2 = FALSE
  if (is.null(access(physeq, "sam_data"))) {
    # Dummy sample data
    SD = sample_data(data.frame(sampnames=sample_names(physeq), sampnames2=sample_names(physeq)))
    rownames(SD) = sample_names(physeq)
    # Add the dummy sample data
    sample_data(physeq) <- SD
    remove.SD = TRUE #use later to remove dummy sample data before function return
  }else if (dim(sample_data(physeq))[2] == 1){
    #If there is a sample_data slot but with only one column, create an additional dummy column - function merge_phyloseq() downstream requires
    sample_data(physeq)[,2] = sample_names(physeq)#dummy column with sample names
    remove.SD2 = TRUE
  }
  
  #Check for taxonomy table, fail if none
  if (is.null(access(physeq, "tax_table"))) {
    stop("The tax_glom.kv() function requires that physeq contain a taxonomyTable")
  }
  #check if there is a phylogenetic tree in the phyloseq object - it needs to be removed first otherwise the function will break trying to 'merge' the tree
  if (!is.null(access(physeq, "phy_tree"))) {
    print("Removing phylogenetic tree")
    physeq@phy_tree <- c()
  }
  tax.k <- data.frame(tax_table(physeq))#extract tax table from phyloseq obj
  #check otu table orientation, change if necessary so that taxa are rows
  if (length(intersect(colnames(otus),rownames(tax.k))>0)){#OTU/ASVs should be rows in otu table, not columns
    otus <-  t(otus)
  }
  tax.k <- tax.k[c(rownames(otus)),]#reorder to match OTUs in specified otu table
  
  levels = colnames(tax_table(physeq))#available taxonomic levels to consider
  for(i in 1:dim(tax.k)[1]){
    x = length(levels)#x=7 (Species column)
    tax.k$lowest[i] <- ifelse(tax.k[i,x]==""|tax.k[i,x]=="NA"|is.na(tax.k[i,x])==TRUE,
                              ifelse(tax.k[i,x-1]==""|tax.k[i,x-1]=="NA"|is.na(tax.k[i,x-1])==TRUE,
                                     ifelse(tax.k[i,x-2]==""|tax.k[i,x-2]=="NA"|is.na(tax.k[i,x-2])==TRUE,
                                            ifelse(tax.k[i,x-3]==""|tax.k[i,x-3]=="NA"|is.na(tax.k[i,x-3])==TRUE,
                                                   ifelse(tax.k[i,x-4]==""|tax.k[i,x-4]=="NA"|is.na(tax.k[i,x-4])==TRUE, "Phylum",
                                                          "Class"),#class
                                                   "Order"),#order
                                            "Family"),#family
                                     "Genus"),#genus
                              "Species")#species
    #so now tax.k should be a vector of the lowest available level of annotation for each OTU (e.g. "Phylum", "Species, "Species", "Class"...)
  }
  #create empty list to store phyloseq sub objects (list of lists)
  phy.list <- list()
  #get the levels that are actually present in tax.k$lowest
  l <- unique(tax.k$lowest)
  for(i in 1:length(l)){
    tax.k.level <- tax.k[tax.k$lowest==l[i],]#subset tax.k by a specific taxonomic rank
    phy.k.level <- prune_taxa(rownames(tax.k.level),physeq)#now subset the phyloseq object to contain only the otus in tax.k.level
    phy.k.merged <- tax_glom(phy.k.level,taxrank=l[i])#now merge at levels[i]
    phy.list[i] <- phy.k.merged
  }
  #now take the list of lists and merge it into one phyloseq object
  i=1
  while(i < length(phy.list)){
    if(i==1){
      phy.new <- merge_phyloseq(phy.list[[i]], phy.list[[i+1]])#first round of merging to create new object named phy.new
    }else {
      phy.new <- merge_phyloseq(phy.new, phy.list[[i+1]])#then subsequently append to phy.new
    }
    i = i+1
  }
  #remove dummy sample data where applicable
  if (remove.SD==TRUE) {
    sample_data(phy.new) <- c()
  }
  if (remove.SD2==TRUE) {
    sample_data(phy.new) <- sample_data(phy.new)[,1]
  }
  #
  print(paste("There are now",ntaxa(phy.new),"merged taxa"))
  return(phy.new)
}

#**********************************
#super.fitZig.kv(): FUNCTION BUILT AROUND METAGENOMESEQ'S FITZIG() AND MRFULLTABLE() FUNCTIONS INCLUDING HEATMAP OF SIGNIFICANT RESULTS
#NB - CURRENTLY ONLY SETUP TO HANDLE TWO-CLASS CATEGORICAL COMPARISONS
#**********************************
#----------------------------------
#ARGUMENTS:
#physeq = phyloseq object with raw counts but pruned by read count (e.g exclude samples with < 10 000 reads)
#factor = the factor of interest (This should be a column name in sample_data(physeq))
#FileName = specify comparison made, not full file path (e.g. "M.prune_PICRUSt_groups_C_vs_A")
#p = adjusted p-value threshold to use when filtering significant results.
#default p=0.05
#perc = threshold for filtering by fraction of +ve ssamples for a given OTU in either group (e.g perc = 0.3 will keep only OTUs where at least one of the two groups have ? 30% of samples +ve for that OTU) 
#default perc = 0.5
#FC = absolute coefficient threshold to use for filtering (default = 1.5)
#main = title of heatmap
#subt = subtitle for heatmap (most commonly filtering settings e.g."FDR < 0.05,|coeff| >= 1.5, >30%+ in either group"
#heatmap.describer = heatmap filename segment specifying clustering or not as well as type of annotation (OTU ID/taxonomy) (e.g."category_ordered_tax_annot.pdf" )
#order = should heatmap columns be sorted manually? TRUE/FALSE, default=TRUE	
#rows.sortby = would you like to sort the resulting heatmap by a column from the rsults table e.g. "adjPavalues" or "coeff"
#labrow = should the heatmap be annotated with taxonomic annotations for OTU ids (TRUE/FALSE; TRUE = taxonomic annotations; FALSE = OTU IDs)
#extra.cols = any extra columns that you would like to plot on the heatmap e.g. c("feeding","HIV_exposure")
#cexCol=scaling factor for heatmap column labels (usually sample names)
#cexRow=scaling factor for heatmap row labels (usually taxon names)
#merged: default FALSE; if data has first been merged supply relevant tax level column name e.g. "Genus"
#-----------------------------------
super.fitZig.kv <- function(physeq, factor, covariate=NULL,outDir,FileName,heatmap.descriptor=NULL,
                            main=NULL, subt=NULL, ordered=TRUE, rows.sortby = NULL,p=0.05, FC = 1.5, perc=0.3, merged = FALSE, extra.cols = NULL,colours=NULL, cexCol=1, cexRow=1){
  #-------------------------------
  #...........................
  #remove any samples with no data for the factor (or covariate)
  sub.index <- which(sample_data(physeq)[,factor] !='NA'& sample_data(physeq)[,factor] !='<NA>'& sample_data(physeq)[,factor] !=''& sample_data(physeq)[,factor] !=' ')
  l = dim(sample_data(physeq))[1]
  if(length(covariate)!=0){
    sub.index.2 <- which(sample_data(physeq)[,covariate] !='NA'& sample_data(physeq)[,covariate] !='<NA>'& sample_data(physeq)[,covariate] !=''& sample_data(physeq)[,covariate] !=' ')
    keep = intersect(sub.index, sub.index.2)
    physeq <- prune_samples(sample_names(physeq)[keep],physeq)
    print(paste(l-length(keep),"of",l,"samples were removed due to missing data"))
  }else{
    physeq <- prune_samples(sample_names(physeq)[sub.index],physeq)
    print(paste(l-length(sub.index),"of",l,"samples were removed due to missing data"))
  }
  #...........................
  #convert phyloseq object to metagenomeSeq object
  MGS <- make_metagenomeSeq(physeq)#
  #...........................
  #Check variables
  levels = levels(pData(MGS)[,factor])
  if(length(covariate)==1){
    if(is.character(pData(MGS)[,covariate])==TRUE){#error message if covariate is not either numeric OR a factor variable
      stop(print("the covariate",covariate,"is a character variable, please convert to factor or numeric first"))
    }		
  }
  if(length(levels)==2){#if we have a factor variable with two levels 
    print(paste(factor,"will be modeled as a binary categorical predictor variable"))
    f = TRUE
  }else if(length(levels) > 2){
    stop("invalid: your factor has more than two levels")#exit function if factor specified does not have exactly two levels
  }else if(length(levels)==0){#then not a factor variable
    if(is.numeric(pData(MGS)[,factor])==TRUE){
      print(paste(factor,"will be modeled as a continuous predictor variable"))
      f=FALSE
    }else if(is.character(pData(MGS)[,factor])==TRUE){#then need to know if data should be numeric or categorical
      stop(paste(factor,"is a character variable, please convert to numeric or factor data first"))
    }
  }
  #DATA NORMALISATION using the cumNormStat function from the metagenomeSeq package, default settings
  par(mar=c(1,1,1,1))
  MGSp = cumNormStat(MGS,pFlag=TRUE,main="Trimmed data")
  MGS = cumNorm(MGS,p=MGSp)
  #-------------------------------
  #-------------------------------
  #Build model based on info specified in the variable 'factor'
  F=pData(MGS)[,factor]
  if(length(covariate)!=0){
    cov.=pData(MGS)[,covariate]
    mod = model.matrix(~F+cov.)#including covariate
  }else{mod = model.matrix(~F)}
  rownames(mod) <- rownames(pData(MGS))
  #-------------------------------
  #fitZig 
  settings = zigControl(maxit=100,verbose=TRUE)
  fit = fitZig(obj = MGS,mod=mod,
               control=settings) #by default, the normalising factors for obj=MGS are included in the model
  #-------------------------------
  #CATEGORICAL VARIABLE
  if(f==TRUE){
    stats <- MRfulltable(fit, coef=colnames(mod)[2],by=colnames(mod)[2],number = dim(fData(MGS))[1],group=2)#data reported for second column in the model
    stats = stats[!is.na(rownames(stats)), ]# if any OTUs left out, rm those from x. Detected by NA rownames.	
  }else if(f==FALSE){
    #OR CONTINUOUS VARIABLE
    stats <- MRcoefs(fit, coef=colnames(mod)[2],by=colnames(mod)[2], number = dim(fData(MGS))[1], group=2)#		
  }
  #-------------------------------	
  
  #NEXT ADD TAXONOMIC DATA FOR SIGNIFICANT OTUs
  #check for NA taxonomic columns, e.g. if merged at genus level, species column will be NAs
  na_counts <- apply(tax_table(physeq),2,function(x) length(which(is.na(x))) == dim(tax_table(physeq))[1])
  if(length(which(na_counts==TRUE))>0){#If there is an NA column(s), remove these
    na_cols <- which(na_counts==TRUE)
    tax_table(physeq) <- tax_table(physeq)[,-c(na_cols)]
  }
  tax.tab <- data.frame(tax_table(physeq))
  tax <- tax.tab[c(rownames(stats)),]
  stats <- cbind(stats, tax)
  #change 'NA' column to coeff.
  n <- which(is.na(colnames(stats)))
  colnames(stats)[n] <- 'coeff'
  #FILTER RESULTS:
  #A) p-values filter
  if(f==TRUE){
    stats.sig <- stats[stats$adjPvalues <= p| stats$fisherAdjP <= p,]
  }else if(f==FALSE){
    stats.sig <- stats[stats$adjPvalues <= p,]
  }
  #B) FC filter
  stats.sig <- stats.sig[abs(stats.sig$coeff) >= FC,]
  nas <- grep("NA",rownames(stats.sig))
  if(length(nas)>0){
    stats.sig = stats.sig[-c(nas), ]# if any OTUs left out, rm those from x. Detected by NA rownames.
  }
  #print results summary to screen:
  if(f==FALSE){
    cat(paste("There were ",dim(stats.sig)[1],"OTUs/ASVs significantly different by",factor, 
              "that met",'\n',"threshold criteria of p",p,"absolute FC",FC))
  }
  
  #GROUP SUMMARIES FOR CATEGORICAL VARIABLES
  if(f==TRUE){
    #only keep results where at least one of the 2 groups has a threshold proportion of + samples:
    N1 <- length(which(mod[,2]==1))#
    N2 <- length(which(mod[,2]==0))#
    N1+N2#
    g1 <- grep("group 1",colnames(stats.sig)[1:2])#
    g0 <- grep("group 0",colnames(stats.sig)[1:2])#
    
    #FILTER BY SPECIFIED % PRESENCE IN AT LEAST ONE GROUP (DEFAULT=30%)
    stats.sig<- stats.sig[stats.sig[,g1] >= perc*N1 | stats.sig$fisherAdjP <= p | stats.sig[,g0] >= perc*N2,]
    nas <- grep("NA",rownames(stats.sig))
    if(length(nas)>0){
      stats.sig = stats.sig[-c(nas), ]# if any OTUs left out, rm those from x. Detected by NA rownames.
    }
    dim(stats.sig)
    #CHANGE "+samples in group1/0" columns to instead reflect mean count calculated across positive samples.
    g1.1 <- grep("group 1",colnames(stats.sig)[3:4])
    g1.1 <- g1.1+2
    g0.1 <- grep("group 0",colnames(stats.sig)[3:4])
    g0.1 <- g0.1+2
    for(i in 1:dim(stats.sig)[1]){
      stats.sig [i,g0.1] <- round(stats.sig[i,g0.1]/stats.sig[i,g0],0)
    }
    colnames(stats.sig)[g0.1] <- "mean_positive_group0"
    for(i in 1:dim(stats.sig)[1]){
      stats.sig [i,g1.1] <- round(stats.sig[i,g1.1]/stats.sig[i,g1],0)
    }
    colnames(stats.sig)[g1.1] <- "mean_positive_group1"	
    stats.sig$percent_positive_group1 <- round(stats.sig[,g1]/N1*100,1)
    stats.sig$percent_positive_group0 <- round(stats.sig[,g0]/N2*100,1)
    #shift these columns to begining of dataframe
    stats.sig <- stats.sig[,c(ncol(stats.sig),(ncol(stats.sig)-1),3:ncol(stats.sig)-2)]
    cat(paste("There were ",dim(stats.sig)[1],"OTUs/ASVs significantly different between",levels[[1]],"vs.",levels[[2]], 
              "that met",'\n',"threshold criteria of p",p,"absolute FC",FC,"and percentage presence in at least one group of",100*perc,"% \n"))
  }
  #-------------------
  #write to file
  #SELECT FILENAME
  #---------------------------------------------------------
  #WRITE RESULTS TO FILE
  file = paste0(outDir,"/",FileName,".csv")	
  #change 'NA' column to coeff.
  print("writing results and model to file")
  write.table(stats.sig,file, sep =",", col.names = NA)
  save(stats.sig, file=paste0(outDir,"/",FileName,".RData"))
  #append model
  finalMod=fit@fit$design
  suppressWarnings(write.table(finalMod,file,append=TRUE , sep =",", col.names = NA))#append model to  file
  #---------------------------------
  #GENERATE HEATMAPS
  #---------------------------------
  #check if any otus to plot
  if(dim(stats.sig)[1]<=1){
    stop("At least two signficant results are required for heatmap plotting")
  }
  #parameter specifications
  file = paste0(outDir,"/",FileName,"_",heatmap.descriptor,".pdf")
  print(file)
  #annotation columns to plot on top of heatmap
  cols = c(factor,extra.cols)
  #first standardize physeq for use in heatmap
  total = median(sample_sums(physeq))
  standf = function(x, t=total) round(t * (x / sum(x)))
  p.std = transform_sample_counts(physeq, standf)
  #manually sort columns by group?
  if(ordered==TRUE){
    heatmap.k(physeq= p.std, data.subset  = stats.sig, subset.samples = rownames(mod), #manual sort columns
              annot.cols = cols,main = main, rows.sortby = rows.sortby, subt=subt,filename = file, order.by = factor, Colv = NA,
              colours=colours, labrow = TRUE, cexCol=cexCol, cexRow=cexRow, merged = merged)
  }else{heatmap.k(physeq= p.std, data.subset  = stats.sig, subset.samples = rownames(mod), #cluster columns
                  annot.cols = cols,main = main, rows.sortby = rows.sortby,subt=subt,filename = file,
                  colours=colours, labrow = TRUE,cexCol=cexCol, cexRow=cexRow, merged = merged)}
  print("making heatmap of results")
  suppressWarnings(warning("super.fitZig.kv"))
  return(stats.sig)
}

#**********************************
#corr.k(): FUNCTION TO PERFORM PAIRWISE CORRELATIONS, TEST FOR SIGNIFICANCE AND PLOT CORRELLOGRAM
#**********************************
#ARGUMENTS:
#mat1: matrix 1 (columns of mat1 to be correlated with columns of mat2)
#mat2: matrix 2
#p: adjusted p-value to sue as cutoff
#use: pairwise correlations? If yes: "pairwise" if not 
#pdf.height: specify height of correlogram .pdf in inches, default=7
#pdf.width: specify width of correlogram .pdf in inches, default=7
corr.k <- function(mat1, mat2, p=0.05, use = "pairwise", adjust = "BH", r = 0.3, filename, rectangles = NULL, pdf.height=NULL, pdf.width=NULL){
  cor.mat = cor(mat1, mat2)
  cor.p <- corr.test(mat1,mat2,use = use,adjust=adjust)#MTC adjustment with corr.test from package 'psych'
  #NOW GET ONLY THE COLUMNS WHERE THERE ARE SIGNIFICANT ENTRIES
  p_mat <- cor.p$p
  if(identical(mat1,mat2)==TRUE){
  diag(cor.mat) <- NA#if mat1=mat2 then we need to exclude diagonals
  diag(p_mat) <- NA#if mat1=mat2 then we need to exclude diagonals
  }
  sig.mat <- cor.mat[,which(apply(p_mat,2,min,na.rm=TRUE) <= p)]
  sig.mat <- sig.mat[which(apply(p_mat,1,min,na.rm=TRUE) <= p),]
  print(paste("There were",dim(sig.mat)[2],"significant after MTC"))
  #also filter on the magnitude of correlation (? 0.3)
  sig.mat <- sig.mat[,which(apply(abs(sig.mat),2,max,na.rm=TRUE) >= r)]
  sig.mat <- sig.mat[which(apply(abs(sig.mat),1,max,na.rm=TRUE) >= r),]
  print(paste("There were",dim(sig.mat)[2],"significant entries after MTC with R >= ",r))
  #select corresponding p-values
  sig.p <- cor.p$p[rownames(sig.mat),colnames(sig.mat)]
  if(dim(sig.mat)[2]<=1){
    return(sig.mat)
  }
  #------------------------------
  #GET CONFIDENCE INTERVALS 
  CIs = na.omit(cor.p$ci[abs(cor.p$ci$r) >= r & abs(cor.p$ci$p) <=p, ])
  #write CIs to file
  write.csv(CIs,file = paste0(outDir,"/",filename,"CIs.csv"))
  #------------------------------
  #now sort by direction of correlation for ease of plot interpretation 
  if(identical(mat1,mat2)==FALSE){
    o.c <- colnames(sig.mat)[order(apply(sig.mat,2,mean), decreasing=FALSE)]
    o.r <- rownames(sig.mat)[order(apply(sig.mat,1,function(x) sum(abs(x))), decreasing=FALSE)]
    sig.mat <- sig.mat[o.r,o.c]
    sig.p <- sig.p[o.r,o.c]
    dim(sig.p)
    #PLOT RESULTS AS CORRELOGRAM
    pdf(paste0(outDir,"/",filename,".pdf"),width=pdf.width, height=pdf.height)#
    corrplot(sig.mat, method='circle',p.mat = sig.p,sig.level=p,
             insig = "blank",tl.cex=0.7,tl.col = "black", cl.cex = 0.7)
    dev.off()
  }else{#PLOT RESULTS AS SQUARE CORRELOGRAM WITH HIERARCHICAL CLUSTERING
    #reset NA diags to 1 for plotting
    diag(sig.mat) <- 1
    pdf(paste0(outDir,"/",filename,".pdf"), width=pdf.width, height=pdf.height)#
    corrplot(sig.mat, method='circle',p.mat = sig.p,sig.level=p,
             insig = "blank",order = "hclust",addrect = rectangles,tl.cex=0.7,tl.col = "black", cl.cex = 0.7)
    dev.off()
  }
}

#**********************************
#MC.summary(): take the output from PICRUSt's metagenome_contributions.py, together with taxonomic annotation for the OTUs included in 
#this table and provides a summary of the contribution of each Family/Genus.. etc to ONE SPECIFIC KEGG gene e.g. K02030
#**********************************
#ARGUMENTS:
#metagenome.contributions = output from PICRUSt's metagenome_contributions.py (imported as a data frame)
#KEGG.gene: KEGG gene of interest e.g. "K02030"
#samples: list of sample names to be included
#level: desired tax level for summary e.g. "Family", "Genus" depending on the column names of your tax table
#outDir: output directory
#file.ext: descriptive file extension such as "Family_K02030"
#tax.table: taxonomic table for OTUs of interest
MC.summary = function(metagenome.contributions, KEGG.gene, tax.table,samples, level,outDir, file.ext){
  out = c()
  for (i in 1:length(samples)) {#for each sample, summarise each family's contribution to the KEGG gene of interest
    Kc = metagenome.contributions[metagenome.contributions$Gene == KEGG.gene,]#subset input by KEGG gene of interest
    this.data=Kc[Kc$Sample==samples[i],]#get data for that specific sample
    otus=this.data[,3]#get all otus for this sample and this gene
    this.tax.table=tax.table[as.character(unique(otus)),]#get corresponding taxonomy for these OTUs.
    tot.taxa=sum(this.data[,7] )#get the percentage contributed by all OTUs in sample i for this specific gene, which should = 1
    taxa=unique(this.tax.table[,level] )#get list of unique taxa at required tax level
    taxa[taxa==" "] = "other"#for taxa with no family level annotation - these will all be summarised together as 'other'
    for (j in 1:length(taxa) ) {
      map2fam = rownames(this.tax.table)[this.tax.table[,level] == taxa[j]]#get OTU IDs for a particular family
      taxa.tot=sum(this.data[this.data$OTU %in% map2fam,7] )/tot.taxa#get the total percent contribution for each family to this specific gene for sample i
      line = t(c(as.character(samples[i]), taxa[j], taxa.tot ))
      out = rbind(line,out)
    }
    print(i)
  }
  #write results to file
  colnames(out) = c("Sample","Taxa","Contribution")
  out = data.frame(out)
  out$Contribution = as.numeric(as.character(out$Contribution))
  write.table(out,file=paste0(outDir,file.ext,'.txt'), quote=FALSE,sep = "\t", col.names=TRUE, row.names=FALSE)
  return(out)
}

#**********************************
#RF.k(): performs random forests analysis on the otu table of a supplied phyloseq object.
#NB: this function is only setup for two-class categorical response variables NOT regression on continuous response variables
#The data is randomly divided into a training (two thirds of the data) and test set (remaining one third of the data not used for training)
#results printed to screen include most important taxa, AUC, PPV, NPV
#Option to specify the top 'x' taxa to see how they perform
#**********************************
#ARGUMENTS:
#data = phyloseq object of interest
#var = variable of interest sample_data(data) column name
#outDir = output directory for results
#train.control: TRUE or FALSE, should RF algorithm parameters be tuned by means of cross-validation (THIS INCREASE RUN TIME, DEFAULT=FALSE)
#cv.fold = cross-validation fold for tuning the RF model (default=10)
#cv.repeat = For repeated k-fold cross-validation only: the number of complete sets of folds to compute
#Nfeatures.validation: 'x' number of top taxa to test (e.g. how good are the top 3 most important taxa at classifying)
#outDir = output directory for results
#testAndtrain: should the dataset be divided into training (randomly select two thirds of the data) and test sets (the remaining one third of the data)
#testAndtrain valid options: TRUE | FALSE
#descriptor: any additional description that needs to go in the file name (e.g. 'merged_otus')
#np: number of top features to display in table and variable importance plot
#positive.class: the positive class needs to be specified - this should match one of the levels() under 'var' of interest e.g:
#for the variable "TB status" the positive class might be "TBpos" - this is used to calculate PPV, NPV
#train.fraction: what fraction of samples should be used for training if testAndtrain=TRUE (the remainder will be used for testing)

RF.k <- function(data, seed = 2,var, cv.fold=10, train.control=FALSE,cv.repeat=3,outDir,Nfeatures.validation=NULL, testAndtrain=FALSE,np=50,
                 positive.class,descriptor = c(),train.fraction=2/3,...){#function default values supplied, change as required
  
  summary_file <- paste0(outDir,"/RF",var,"_results_",cv.fold,"_", Nfeatures.validation,"_",seed,descriptor,".txt")
  write(print(paste(length(which(is.na(sample_data(data)[,var]))),"samples did not have response variable data, removing these...")), 
        file = summary_file,sep = "\t")
  sub.index <- which(!is.na(sample_data(data)[,var]))#
  data <- prune_samples(sample_names(data)[sub.index],data)
  if(testAndtrain==TRUE){
    table = otu_table(data)
    train.n = round(train.fraction*nsamples(data))
    #generate random integers between 1 and nsamples(phy.temp)
    set.seed(seed)
    train.rows=createDataPartition(y=unlist(sample_data(data)[,var]), p=train.fraction, list=FALSE)#randomly select desired fraction of samples for training and testing while balancing class distributions within splits
    train = table[,train.rows]
    test = table[,-train.rows]
    tbl <- table(sample_data(data)[colnames(train),var])
    write(print(paste("Training set size: ",dim(train)[2], "samples","with",tbl[1],"and",tbl[2], "samples per class")), 
          file = summary_file,sep = "\t", append=T)
    tbl <- table(sample_data(data)[colnames(test),var])
    write(print(paste("Test set size: ",dim(test)[2], "samples","with",tbl[1],"and",tbl[2], "samples per class")),
          file = summary_file,sep = "\t", append=T)
    #TRAINING - RANDOMLY SELECT 'N=training.fraction' OF DATASET
    # Make a dataframe of training data with OTUs as column and samples as rows
    predictors <- t(train)
    predictors = data.frame(predictors)
    #------------------
    # Make one column for our outcome/response variable 
    response <- as.factor(unlist(sample_data(data)[train.rows,var]))
    names(response) <- sample_names(data)[train.rows]
  }
  else if(testAndtrain==FALSE){
    table = otu_table(data)
    tbl <- table(sample_data(data)[,var])
    write(print(paste("Data set size: ",tbl[1]+tbl[2], "samples","with",tbl[1],"and",tbl[2], "samples per class")), 
          file = summary_file,sep = "\t", append=T)
    # Make a dataframe of training data with OTUs as column and samples as rows
    predictors <- t(table)
    predictors = data.frame(predictors)
    #------------------
    # Make one column for our outcome/response variable 
    response <- as.factor(unlist(sample_data(data)[,var]))
    names(response) <- sample_names(data)
  }
  
  #------------------------------------
  #BUILD RF MODEL
  # Combine them into 1 data frame
  rf.data <- data.frame(response, predictors)
  head(rf.data)
  set.seed(seed)
  
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(print("******************************************************************************"),file = summary_file,sep = "\t", append=T)
  write(print("Building RF model on training set....."),file = summary_file,sep = "\t", append=T)
  if(train.control==TRUE){
    write(print(paste("Performing model tuning:",cv.repeat,"round(s) of",cv.fold,"fold cross-validation")),file = summary_file,sep = "\t", append=T)
    write(print("******************************************************************************"),file = summary_file,sep = "\t", append=T)
    
    classify<-train(response~.,data=rf.data,method="rf",
                    trControl=trainControl(method="repeatedcv",repeats=cv.repeat,number=cv.fold), metric="Accuracy", tuneLength=10)
    print(classify)
  }
  else if(train.control==FALSE){
    write(print("Performing model tuning.."),file = summary_file,sep = "\t", append=T)
    write(print("******************************************************************************"),file = summary_file,sep = "\t", append=T)
    classify<-train(response~.,data=rf.data,method="rf", metric="Accuracy")
    print(classify)
  }
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Cmat = confusionMatrix(data = classify$finalModel$predicted, reference = rf.data$response, positive=positive.class)
  print(Cmat)
  sink(summary_file, append=T)
  print(Cmat)
  sink()
  #------------------------------------------
  print("******************************************************************************")
  write(print("******************************************************************************"),file = summary_file,sep = "\t", append=T)
  #------------------------------------------
  #RF CV for feature selection:
  rf.cv = rfcv(rf.data[,2:dim(rf.data)[2]], rf.data[,1], cv.fold=cv.fold)
  write(print("Cross-validated error rates associated with stepwise reduction of features:"),file = summary_file,sep = "\t", append=T)
  print(rf.cv$error.cv)
  sink(summary_file, append=T)
  print(rf.cv$error.cv)
  sink()
  pdf(paste0(outDir,"/RF_elbow_plot_",var,"_",cv.fold,"_",Nfeatures.validation,"_",seed,descriptor,".pdf"))
  with(rf.cv, plot(n.var, error.cv, log="x", type="o", lwd=2))
  dev.off()
  #----------------------------------
  # Make a data frame with predictor names and their importance
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  imp=importance(classify$finalModel) #When we use metric="Accuracy" without importance=TRUE in the model build train() we get the same results we would've when using the randomForest package (meanDecreaseGini)
  imp <- data.frame(predictors = rownames(imp), imp)
  # Order the predictor levels by importance
  imp.sort <- arrange(imp, desc(MeanDecreaseGini))
  imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
  # Select the top n predictors
  imp.n <- imp.sort[1:np, ]
  #add taxonomy info
  rownames(imp.n) = imp.n[,1]
  labs = tax.lab(imp.n, data, labrow=TRUE,merged=FALSE)
  imp.n$tax = labs
  imp.n$tax = factor(imp.n$tax, levels = imp.n$tax)#avoid ggplot sorting alphabetically\
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(print("******************************************************************************"),file = summary_file,sep = "\t", append=T)
  write(print(paste("THE TOP",np," MOST IMPORTANT FEATURES WERE:")),file = summary_file,sep = "\t", append=T)
  print(imp.n)
  sink(summary_file, append=T)
  print(imp.n)
  sink()
  write(print("******************************************************************************"),file = summary_file,sep = "\t", append=T)
  p <- ggplot(imp.n,aes(x = imp.n[,"tax"], y = imp.n[,"MeanDecreaseGini"]),environment = environment())+
    geom_bar(stat = "identity", fill = "grey") +
    coord_flip() +
    ggtitle(paste("Most important taxa for classifying samples by ",var))+
    theme_bw() + 
    theme(panel.background = element_blank(), 
          panel.grid.major = element_blank(),  #remove major-grid labels
          panel.grid.minor = element_blank(),  #remove minor-grid labels
          plot.background = element_blank())
  pdf(paste0(outDir, "/RFs",var,"variable_importance_plot_",cv.fold,"_",Nfeatures.validation,"_",seed,descriptor,".pdf"))#note: ggplot doesn't plot within a function - this is a workaround using print()
  print(p)
  dev.off()
  #------------------------------------
  
  #SEE HOW THE TOP X NUMBER OF FEATURES DO - SPECIFY IN 'Nfeatures.validation' PARAMETER
  if(length(Nfeatures.validation)!=0){
    #sort by mean decrease in Gini index and select top 'x'
    goodPredictors = imp.n[1:Nfeatures.validation,1]
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    classify.x<-train(response~.,data=rf.data[,c(as.character(goodPredictors),"response")],method="rf",metric="Accuracy")
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(print("******************************************************************************"),file = summary_file,sep = "\t", append=T)
    write(print(paste("TRAINING SET classification summary if using the TOP",Nfeatures.validation, "features only")),summary_file,sep = "\t", append=T)
    write(print("******************************************************************************"),file = summary_file,sep = "\t", append=T)
    write(print(paste("Feature(s) selected:",goodPredictors)),file = summary_file,sep = "\t", append=T)
    print(classify.x$finalModel)
    #confusion matrix:
    Cmat.x = confusionMatrix(data = classify.x$finalModel$predicted, reference = rf.data$response, positive=positive.class)
    print(Cmat.x)
    sink(summary_file, append=T)
    print(Cmat.x)
    sink()
    
  }
  if(testAndtrain==TRUE){
    #---------------------------------------
    #NOW TEST CLASSIFIER IN TEST COHORT (I.E. REMAINING 1/3 OF DATASET) USING THE TOP X FEATURES
    #---------------------------------------
    test <- t(test)
    test = data.frame(test)
    
    if(length(Nfeatures.validation)!=0){
      rf.test = classify.x
    }
    else{rf.test = classify}
    #Predict classes based on model built on training data (reduced features set)
    predict.test = predict(rf.test, test)
    #Get actual classes for response variable, test samples
    test.response <- as.factor(unlist(sample_data(data)[-train.rows,var]))
    names(test.response) <- sample_names(data)[-train.rows]
    #identical(names(test.response),rownames(test))
    Cmat.test=confusionMatrix(data=predict.test, reference=test.response, positive=positive.class)
    
    write(print("******************************************************************************"),file = summary_file,sep = "\t", append=T)
    write(print(paste("TEST SET model testing (",round((1-train.fraction)*100) ,"% of the data )")),file = summary_file,sep = "\t", append=T)
    if(length(Nfeatures.validation)!=0){
      write(print(paste("using the top",Nfeatures.validation,"features")),file = summary_file,sep = "\t", append=T)
    }
    else{write(print("using all features"),file = summary_file,sep = "\t", append=T)
    }
    write(print("******************************************************************************"),file = summary_file,sep = "\t", append=T)
    print(Cmat.test)
    sink(summary_file, append=T)
    print(Cmat.test)
    sink()
    #---------------------------------
    
  }
}


#THE END
#///////////////////////////////////////

# > sessionInfo()
#R version 3.6.0 (2019-04-26)
#Platform: x86_64-apple-darwin15.6.0 (64-bit)
#Running under: macOS Sierra 10.12.6

#Matrix products: default
#BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
#LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

#Random number generation:
# RNG:     Mersenne-Twister 
# Normal:  Inversion 
# Sample:  Rounding 

#locale:
#[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# attached base packages:
# [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:

#other attached packages:
# [1] gridExtra_2.3        caret_6.0-84         ROCR_1.0-7           gplots_3.0.1.1       dplyr_0.8.1          randomForest_4.6-14  metagenomeSeq_1.26.3
# [8] RColorBrewer_1.1-2   glmnet_2.0-18        foreach_1.4.4        Matrix_1.2-17        limma_3.40.2         matrixStats_0.54.0   psych_1.8.12        
#[15] corrplot_0.84        ggplot2_3.2.0        vegan_2.5-5          lattice_0.20-38      permute_0.9-5        NMF_0.23.6           Biobase_2.44.0      
#[22] BiocGenerics_0.30.0  cluster_2.0.9        rngtools_1.3.1.1     pkgmaker_0.28        registry_0.5-1       phyloseq_1.28.0       
#///
