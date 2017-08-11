library(phyloseq)
library(NMF)
library(vegan)
library(ggplot2)
library(matrixStats)#rowSds
library(fifer)
library(metagenomeSeq)#differential abundance testing
library(randomForest)
library(dplyr)
library(ROCR)
library(gridExtra)#for grid.arrange()
nmf.options(grid.patch=TRUE)#set to avoid blank first pdf page being created

#DEFINE CUSTOM COLOR PALETTE WHERE COLORS ARE EASIER TO DISTINGUISH FROM ONE ANOTHER
myPalette <- c('#89C5DA', "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861")

#--------------------------------------
#THIS METAGENOME SEQ FUNCTION WAS COPIED AND MODIFIED FROM CODE USED FOR THE MCMURDIE 2014 PAPER - http://joey711.github.io/waste-not-supplemental/simulation-differential-abundance/simulation-differential-abundance-server.html
#---------------------------------------
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
#TAX.LAB - GET LOWEST LEVEL OF TAXONOMIC ANNOTATION AVAILABLE TO SUBSTITUTE OTUS FOR PLOTTING PURPOSES.
tax.lab <- function(otus, physeq, labrow=TRUE,merged=FALSE){
	tax <- data.frame(tax_table(physeq))
	tax <- tax[c(rownames(otus)),]
	tax[is.na(tax)] <- ""
	if(merged==FALSE){
		tax$unique <- rep("",dim(tax)[1])#make unique names
		for(i in 1:dim(tax)[1]){
			x = 7#column 7 = Species
			tax$unique[i] <- ifelse(tax[i,x]==""|is.na(tax[i,x])==TRUE,
					ifelse(tax[i,x-1]==""|is.na(tax[i,x-1])==TRUE,
							ifelse(tax[i,x-2]==""|is.na(tax[i,x-2])==TRUE,
									ifelse(tax[i,x-3]==""|is.na(tax[i,x-3])==TRUE,#phylum
											ifelse(tax[i,x-4]==""|is.na(tax[i,x-4])==TRUE, as.character(tax[i,x-5]),
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
			if(is.na(tax$Species[i])==FALSE & tax$Species[i]!=""){#get OTUs with species-level annotation
				get = substr(tax$Genus[i],1,1)
				tax$unique[i] = sub(tax$unique[i],paste0(get,".",tax$unique[i]),tax$unique[i])
			}	
		}
		tax$unique <- make.unique(tax$unique)
#---------------------
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



#otus = otu table of interest (often subset)
#physeq = phyloseq object of interst (for tax annotation)
#labrow: TRUE=annotate rows with taxonomy info, FALSE=annotate rows with OTU info.
#if OTUs have been merged, TRUE, else FALSE
#---------------------
#Generic heatmap function using package NMF
#---------------------
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

#---------------------
#Generic heatmap function using package NMF
#---------------------
heatmap.k = function(physeq=NULL, plot.otus = TRUE,data.subset=NULL,subset.samples = NULL, 
		annot.cols = NULL,colours = NULL, order.by = NULL,main=NULL, subt=NULL,filename,
		Colv = NULL,rows.sortby=NULL,labrow = FALSE, merged = FALSE, distfun = 'euclidean',hclustfun = 'average', cexCol=0.5, cexRow=0.3, scale = "none"){
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
			otus[zs] <- 0.1 #set zs to 0.1
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
	#a.track = a.track[colnames(plot),]
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
		labs <- tax.lab(physeq=physeq,otus=plot,labrow, merged)
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

#---------------------
#Generic barplot function build on phyloseq plot_bar(): NB: THE NUMBER OF SAMPLES IN EACH GROUP IS ONLY DISPLAYED CORRECTLY FOR TWO GROUPS (NOT MORE)
#---------------------
#This function was modified from the phyloseq plot_bar() function where ggplot2's geom_bar no longers sorted stacked bars by abundance.
plot_bar_fix <- function (physeq, x = "Sample", y = "Abundance", fill = level, 
		title = NULL, facet_grid = NULL) 
{
	mdf = psmelt(physeq)
	p=ggplot(mdf[order(-mdf[,y]),],aes_string(x= x, y=y, fill = fill)) + geom_bar(stat='identity')
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
		level, perc, order.bars = NULL,x.axis=NULL, y.axis=NULL, filen = "",outDir){
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
	merged = filter_taxa(merged, function(x) sum(x > count) > (perc*length(x)), TRUE)#filter taxa for barplot display
	new <- merge_samples(merged,cat, fun = "median")#Merge by category
	new <- transform_sample_counts(new, function(x) 100 * x/sum(x))#convert to % abundance
	sample_data(new)[,cat] <- rownames(sample_data(new)) #merge_samples also "merges the sample table, replace relevant cat column with rownames
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
	p2 <- plot_bar_fix(new,x = cat, fill=level ,title = paste0(level,"-level 16S ",x.axis,". Ns = ",for.title))+ coord_flip() + 
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


#tax_glom.kv - combination of phyloseq's tax_glom() and my tax.lab() function
#GET LOWEST LEVEL OF TAXONOMIC ANNOTATION AVAILABLE AND MERGE COUNTS AT THIS LEVEL (i.e. I have 3 otus of L.iners and another 3 of Lactobacillus (unknown species) -->
#e.g. merge L.iners at species level and Lactobacillus (no species annot) at genus level.
#modified from phyloseq function tax_glom
#NB: this works for when we have the standard 7 ranks
#colnames(tax_table(physeq))
#[1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"
#otus = otu table of interest, if not specified, defaults to OTU table of the phyloseq obj. specified.
#physeq = phyloseq object with taxonomic annotation
tax_glom.kv <- function (physeq) 
{
	otus = otu_table(physeq)
	if (is.null(access(physeq, "tax_table"))) {
		stop("The tax_glom.kv() function requires that physeq contain a taxonomyTable")
	}
	#check if there is a phylogenetic tree in the phyloseq object - it needs to be removed first otherwise the function will break trying to 'merge' the tree
	if (!is.null(access(physeq, "phy_tree"))) {
		print("Removing phylogenetic tree")
		physeq@phy_tree <- c()
	}
	tax.k <- data.frame(tax_table(physeq))#extract tax table from phyloseq obj
	tax.k <- tax.k[c(rownames(otus)),]#reorder to natch OTUs in specified otu table
	tax.k[is.na(tax.k)] <- ""#set NAs to ""
	levels = colnames(tax_table(physeq))#available taxonomic levels to consider
	for(i in 1:dim(tax.k)[1]){
		x = length(levels)#x=7 (Species column)
		tax.k$lowest[i] <- ifelse(tax.k[i,x]==""|is.na(tax.k[i,x])==TRUE,
				ifelse(tax.k[i,x-1]==""|is.na(tax.k[i,x-1])==TRUE,
						ifelse(tax.k[i,x-2]==""|is.na(tax.k[i,x-2])==TRUE,
								ifelse(tax.k[i,x-3]==""|is.na(tax.k[i,x-3])==TRUE,
										ifelse(tax.k[i,x-4]==""|is.na(tax.k[i,x-4])==TRUE, "Phylum",
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
	print(paste("There are now",ntaxa(phy.new),"merged taxa"))
	return(phy.new)
}

#------------------------------------------------------
#FUNCTION BUILT AROUND METAGENOMESEQ'S FITZIG() AND MRFULLTABLE() FUNCTIONS INCLUDING HEATMAP OF SIGNIFICANT RESULTS
#NB - CURRENTLY ONLY SETUP TO HANDLE TWO-CLASS CATEGORICAL COMPARISONS
#----------------------------------
super.fitZig.kv <- function(physeq, factor, covariate=NULL,outDir,FileName,heatmap.descriptor=NULL,
		main=NULL, subt=NULL, ordered=TRUE, rows.sortby = NULL,p=0.05, FC = 1.5, perc=0.3, extra.cols = NULL,colours=NULL){
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
		cat(paste("There were ",dim(stats.sig)[1],"OTUs significantly different by",factor, 
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
		cat(paste("There were ",dim(stats.sig)[1],"OTUs significantly different between",levels[[1]],"vs.",levels[[2]], 
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
	finalMod=fit$fit$design
	write.table(finalMod,file,append=TRUE , sep =",", col.names = NA)#append model to  file
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
				colours=colours, labrow = TRUE)
	}else{heatmap.k(physeq= p.std, data.subset  = stats.sig, subset.samples = rownames(mod), #cluster columns
				annot.cols = cols,main = main, rows.sortby = rows.sortby,subt=subt,filename = file,
				colours=colours, labrow = TRUE)}
	print("making heatmap of results")
	suppressWarnings(warning("super.fitZig.kv"))
	return(stats.sig)
}

#------------------------------------------------------
#physeq = phyloseq object with raw counts but pruned by read count (e.g exclude samples with < 10 000 reads)
#factor = the factor of interest (This should be a column name in sample_data(physeq))
#FileName = specify comparison made, not full file path (e.g. "M.prune_PICRUSt_groups_C_vs_A")
#p = adjusted p-value threshold to use when filtering significant results.
#default p=0.05
#perc = threshold for filtering by fraction of +ve ssamples for a given OTU in either group (e.g perc = 0.3 will keep only OTUs where at least one of the two groups have ³ 30% of samples +ve for that OTU) 
#default perc = 0.5
#FC = absolute coefficient threshold to use for filtering (default = 1.5)
#main = title of heatmap
#subt = subtitle for heatmap (most commonly filtering settings e.g."FDR < 0.05,|coeff| >= 1.5, >30%+ in either group"
#heatmap.describer = heatmap filename segment specifying clustering or not as well as type of annotation (OTU ID/taxonomy) (e.g."category_ordered_tax_annot.pdf" )
#order = should heatmap columns be sorted manually? TRUE/FALSE, default=TRUE	
#rows.sortby = would you like to sort the resulting heatmap by a column from the rsults table e.g. "adjPavalues" or "coeff"
#labrow = should the heatmap be annotated with taxonomic annotations for OTU ids (TRUE/FALSE; TRUE = taxonomic annotations; FALSE = OTU IDs)
#extra.cols = any extra columns that you would like to plot on the heatmap e.g. c("feeding","HIV_exposure")
#------------------------------------------------------

#--------------------------------------------------
#The function RF.k performs random forests analysis on the otu table of a supplied phyloseq object.
#NB: this function is only setup for categorical response variables NOT regression on continuous response variables
#The data is randomly divided into a training (two thirds of the data) and test set (remaining one third of the data not used for training)
#results printed to screen include most important taxa, AUC, PPV, NPV
#Option to specify the top 'x' taxa to see how they perform
#--------------------------------------------------
#data = phyloseq object of interest
#var = variable of interest sample_data(data) column name
#outDir = output directory for results
#ntree = number of trees for randomForest() function
#mtry = Number of variables randomly sampled as candidates at each split; default as supplied in the randomForest() function
#cv.fold = cross-fold validation (default=3)
#Nfeatures.validation: 'x' number of top taxa to test (e.g. how good are the top 3 most important taxa at classifying)
#outDir = output directory for results
#testAndtrain: should the dataset be divided into training (randomly select two thirds of the data) and test sets (the remaining one third of the data)
#testAndtrain valid options: TRUE | FALSE
#descriptor: any additional description that needs to go in the file name (e.g. 'merged_otus')
RF.k <- function(data, seed = 212,var, ntree=500, mtry=if (!is.null(y) && !is.factor(y))
					max(floor(ncol(x)/3), 1) else floor(sqrt(ncol(x))), cv.fold=3, outDir,Nfeatures.validation=NULL, testAndtrain=FALSE,descriptor = c(),...){#function default values supplied, change as required
	
	summary_file <- paste0(outDir,"/RF",var,"_results_",ntree,"_",cv.fold,"_", Nfeatures.validation,"_",seed,descriptor,".txt")
	write(print(paste(length(which(is.na(sample_data(data)[,var]))),"samples did not have response variable data, removing these...")), 
			file = summary_file,sep = "\t")
	sub.index <- which(!is.na(sample_data(data)[,var]))#
	data <- prune_samples(sample_names(data)[sub.index],data)
	if(testAndtrain==TRUE){
		table = otu_table(data)
		train.n = round(2/3*nsamples(data))
		#generate random integers between 1 and nsamples(phy.temp)
		set.seed(seed)
		random = sample(1:nsamples(data), train.n)
		train = table[,random]
		test = table[,-random]
		tbl <- table(sample_data(data)[colnames(train),var])
		write(print(paste("Training set size: ",dim(train)[2], "samples","with",tbl[1],"and",tbl[2], "samples per class")), 
				file = summary_file,sep = "\t", append=T)
		tbl <- table(sample_data(data)[colnames(test),var])
		write(print(paste("Test set size: ",dim(test)[2], "samples","with",tbl[1],"and",tbl[2], "samples per class")),
				file = summary_file,sep = "\t", append=T)
		#TRAINING - RANDOMLY SELECT 2/3 OF DATASET
		# Make a dataframe of training data with OTUs as column and samples as rows
		predictors <- t(train)
		predictors = data.frame(predictors)
		#------------------
		# Make one column for our outcome/response variable 
		response <- as.factor(unlist(sample_data(data)[random,var]))
		names(response) <- sample_names(data)[random]
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
	set.seed(2)
	classify <- randomForest(response~ ., data = rf.data, importance=T, proximity=T,ntree = 10000, na.action=na.omit)
	classify
	sink(summary_file,append=T)
	print(classify)
	sink()
	#------------------------------------------
	#RF CV for feature selection:
	rf.cv = rfcv(rf.data[,2:dim(rf.data)[2]], rf.data[,1], cv.fold=cv.fold)
	write(print("Cross-validated error rates associated with stepwise reduction of features:"),file = summary_file,sep = "\t", append=T)
	print(rf.cv$error.cv)
	sink(summary_file, append=T)
	print(rf.cv$error.cv)
	sink()
	pdf(paste0(outDir,"/RF_elbow_plot_",var,"_",ntree,"_",cv.fold,"_",Nfeatures.validation,"_",seed,descriptor,".pdf"))
	with(rf.cv, plot(n.var, error.cv, log="x", type="o", lwd=2))
	dev.off()
	#----------------------------------
	# Make a data frame with predictor names and their importance
	imp <- importance(classify)
	imp <- data.frame(predictors = rownames(imp), imp)
	# Order the predictor levels by importance
	imp.sort <- arrange(imp, desc(MeanDecreaseGini))
	imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
	# Select the top 20 predictors
	imp.20 <- imp.sort[1:20, ]
	#add taxonomy info
	rownames(imp.20) = imp.20[,1]
	labs = tax.lab(imp.20, data, labrow=TRUE,merged=FALSE)
	#head(labs)
	#rownames(imp.20) = labs
	#head(imp.20)
	imp.20$tax = labs
	imp.20$tax = factor(imp.20$tax, levels = imp.20$tax)#avoid ggplot sorting alphabetically\
	write(print("*****************************"),file = summary_file,sep = "\t", append=T)
	write(print("THE TOP 20 MOST IMPORTANT FEATURES WERE:"),file = summary_file,sep = "\t", append=T)
	print(imp.20)
	sink(summary_file, append=T)
	print(imp.20)
	sink()
	write(print("*****************************"),file = summary_file,sep = "\t", append=T)
	p <- ggplot(imp.20,aes(x = imp.20[,"tax"], y = imp.20[,"MeanDecreaseGini"]),environment = environment())+
			geom_bar(stat = "identity", fill = "grey") +
			coord_flip() +
			ggtitle(paste("Most important taxa for classifying samples by ",var))+
			theme_bw() + 
			theme(panel.background = element_blank(), 
					panel.grid.major = element_blank(),  #remove major-grid labels
					panel.grid.minor = element_blank(),  #remove minor-grid labels
					plot.background = element_blank())
	pdf(paste0(outDir, "/RFs",var,"variable_importance_plot_",ntree,"_",cv.fold,"_",Nfeatures.validation,"_",seed,descriptor,".pdf"))#note: ggplot doesn't plot within a function - this is a workaround using print()
	print(p)
	dev.off()
	#------------------------------------
	#------------------------------------
	# Calculate ROC AUC
	predictions = as.vector(classify$votes[,2])
	pred1 = prediction(predictions, classify$y)#
	auc1 = performance(pred1, 'auc')
	#------------------------------------
	#Calculate AUC, PPV, NPV
	write(print(paste0("Training AUC=",round(auc1@y.values[[1]],2))),file = summary_file,sep = "\t", append=T)
	# Calculate predictive value of classifier
	confusion1 = classify$confusion[c(2,1),c(2,1)]
	ppv1 = ((confusion1[1,1])/(confusion1[1,1]+confusion1[2,1]))
	write(print(paste0("Training PPV=",round(ppv1,2))),file = summary_file,sep = "\t", append=T)
	npv1 = ((confusion1[2,2])/(confusion1[1,2]+confusion1[2,2]))
	write(print(paste0("Training NPV=",round(npv1,2))),file = summary_file,sep = "\t", append=T)
	# Plot ROC AUC for classifier 
	pdf(paste0(outDir,"/RFs_",var,"_variable_ROC_curve_training_set_",ntree,"_",cv.fold,"_",Nfeatures.validation,"_",seed,descriptor,".pdf"))
	perf1 = performance(pred1, 'tpr','fpr')
	plot(perf1, main='ROC Curve taxa predictors, training set', col='red', lwd=2)
	text(0.5,0.5,paste('AUC = ',format(auc1@y.values[[1]],digits=5,scientific=FALSE),'\nPPV = ',
					format(ppv1,digits=5,scientific=FALSE),'\nNPV = ',format(npv1,digits=5,scientific=FALSE)),cex=2)
	dev.off()
	#SEE HOW THE TOP X NUMBER OF FEATURES DO - SPECIFY IN 'Nfeatures.validation' PARAMETER
	if(length(Nfeatures.validation)!=0){
		#sort by mean decrease in Gini index and select top 'x'
		goodPredictors = rownames(classify$importance)[order(classify$importance[,4],decreasing=T)][1:Nfeatures.validation]
		rf.x = randomForest(response ~ .,data=rf.data[,c(goodPredictors,"response")],importance=T,proximity=T,ntree = ntree, na.action=na.omit)#~. include all variables
		write(print("*****************************"),file = summary_file,sep = "\t", append=T)
		write(print(paste("Training set classification summary if using the top",Nfeatures.validation, "features only")),summary_file,sep = "\t", append=T)
		write(print(paste("Feature(s) selected:",goodPredictors)),file = summary_file,sep = "\t", append=T)
		print(rf.x$importance[order(rf.x$importance[,4], decreasing=T),])
		sink(summary_file, append=T)
		print(rf.x$importance[order(rf.x$importance[,4], decreasing=T),])
		sink()
		#write(print(paste("Training set classification summary if using the top",Nfeatures.validation, "features only")),file = paste0(outDir,"RF",var,"results.txt"),sep = "\t", append=T)
		print(rf.x)
		sink(summary_file, append=T)
		print(rf.x)
		sink()
	}
	if(testAndtrain==TRUE){
		#---------------------------------------
		#NOW TEST CLASSIFIER IN TEST COHORT (I.E. REMAINING 1/3 OF DATASET) USING THE TOP X FEATURES
		#---------------------------------------
		if(length(Nfeatures.validation)!=0){
			rf.test = rf.x
		}
		else{rf.test = classify}
		predictors <- t(test)#now use test set
		response <- as.factor(unlist(sample_data(data)[-random,var]))
		names(response) <- sample_names(data)[-random]
		test.data <- data.frame(response, predictors)	
		# Getting predictions for testing hold-out data
		validation.resp1 = predict(rf.test, test.data, type='response')
		validation.vote1 = predict(rf.test, test.data, type='vote')
		# Create confusion matrix
		confusion1 = table(data.frame(cbind(Actual=response,Predicted=validation.resp1)))
		class.labels = levels(as.factor(unlist(sample_data(data)[,var])))
		rownames(confusion1) = c(class.labels[1],class.labels[2])
		colnames(confusion1) = c(class.labels[1],class.labels[2])
		ro1 = c(class.labels[1],class.labels[2])
		#predicted error
		err1 = ((confusion1[1,2]+confusion1[2,1])/sum(confusion1))
		write(print("*****************************"),file = summary_file,sep = "\t", append=T)
		write(print("moving to TEST set (1/3 of data)"),file = summary_file,sep = "\t", append=T)
		if(length(Nfeatures.validation)!=0){
			write(print(paste("using the top",Nfeatures.validation,"features")),file = summary_file,sep = "\t", append=T)
		}
		else{write(print("using all features"),file = summary_file,sep = "\t", append=T)
		}
		write(print("*****************************"),file = summary_file,sep = "\t", append=T)
		write(print(paste0("Validation predicted error: ",round(err1*100,2),"%")),file = summary_file,sep = "\t", append=T)
		#--------------------------
		# Calculate ROC AUC for testing set
		validation.predictions = as.vector(validation.vote1[,2])
		validation.pred1 = prediction(validation.predictions, response)
		validation.auc1 = performance(validation.pred1, 'auc')
		write(print(paste0("Test set AUC=",round(auc1@y.values[[1]],2))),summary_file,sep = "\t", append=T)
		ppv1 = ((confusion1[1,1])/(confusion1[1,1]+confusion1[2,1]))
		write(print(paste0("Test set PPV=",round(ppv1,2))),summary_file,sep = "\t", append=T)
		npv1 = ((confusion1[2,2])/(confusion1[1,2]+confusion1[2,2]))
		write(print(paste0("Test set NPV=",round(npv1,2))),summary_file,sep = "\t", append=T)
		#------------------------------------
		pdf(paste0(outDir,'/RFs_',var,'_validation.taxa_',ntree,"_", cv.fold, "_",Nfeatures.validation,"_",seed,descriptor,'_ROC_curve.pdf'))
		validation.perf1 = performance(validation.pred1, 'tpr','fpr')
		plot(validation.perf1, main='ROC Curve taxa predictors, test set', col='blue', lwd=2)
		text(0.5,0.5,paste('AUC = ',format(validation.auc1@y.values[[1]],digits=5,scientific=FALSE),'\nPPV = ',format(ppv1,digits=5,scientific=FALSE),'\nNPV = ',format(npv1,digits=5,scientific=FALSE)),cex=2)
		dev.off()
	}
}

#////////////////////////////////
#///////////////////////////////////////
