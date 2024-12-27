## Metabolite level prediction from gene expression data of Arabidopsis in
## Lozano-Elena et al. 2022 
## by ccle network
## May 24th, 2022                                                                
## *************************************

## Settings ----------------------------
set.seed(1)
rm(list=ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, plyr, readr, tidyverse, fs, DESeq2, corto, openxlsx)	

setwd("~/Documents/metabpredplants/LozanoElena2022")
source("~/Documents/metabpredplants/code/vst.R")

## Cleaning data ----------------------------
# load expression data
files <- dir_ls("~/Documents/metabpredplants/LozanoElena2022/data/GSE119382_RAW", recurse = TRUE, type = "file", glob = "*.txt")  
fname <- "~/Documents/metabpredplants/LozanoElena2022/data/GSE119382_RAW/normcounts.rda"

if(!file.exists(fname)){
	list_exp <- map(files, read_table)
	rawcounts <- (Reduce(function(...) merge(..., by = "gene", all=TRUE), list_exp))
	rownames(rawcounts) <- rawcounts$gene
	rawcounts = rawcounts[,!(colnames(rawcounts) == 'gene')]
	drop_rows <- c('__no_feature', '__ambiguous', '__too_low_aQual',  '__not_aligned', '__alignment_not_unique')
	rawcounts <- subset(rawcounts, !rownames(rawcounts) %in% drop_rows)
	normcounts <- vst(data.matrix(rawcounts))
	save(normcounts, file=fname)
}else{load(fname)} 

## DIOPT orthologs predictions
# load DIOPT output 
load("~/Documents/metabpredplants/LozanoElena2022/data/orthologs/orthologs.rda")

matchids_row <- function(df, oldids, newids) {
	i <- match(oldids,rownames(df), nomatch=0)
	rownames(df)[i] <- as.character(newids)[i]
	return(df)
	}	
		
fname <- "~/Documents/metabpredplants/LozanoElena2022/data/GSE119382_RAW/normcountsids.rda"
if(!file.exists(fname)){
	# replace arabidopsis genes with corresponding orthologs
	normcounts <- subset(normcounts, rownames(normcounts) %in% orthologs$Search.Term)
	
	oldids=orthologs$Search.Term
	newids=orthologs$Human.Symbol
	df=normcounts

	normcountsids <- matchids_row(df, oldids, newids)
	## remove the following exp samples because matching with the same metabolomics sample
	drop_cols <- c("GSM3373994_D_3","GSM3373998_W_1")
	normcounts = normcountsids[, !(colnames(normcountsids) %in% drop_cols)]

	save(normcounts, file=fname)
}else{load(fname)} 
# dim(normcountsids) #[1] 4973    9

# load metabolite data
fname <- "~/Documents/metabpredplants/LozanoElena2022/data/normmetabids.rda"

if(!file.exists(fname)){
	pairedsamples <- read.delim('~/Documents/metabpredplants/LozanoElena2022/data/pairedSamples.txt', sep='\t', header=TRUE)
	metabsamples <- pairedsamples$Metabolomics
	#which(duplicated(metabsamples))
	#[1] 8 
	#15316oA_72
	metab <- read.xlsx("~/Documents/metabpredplants/LozanoElena2022/data/41597_2022_1161_MOESM3_ESM_metab.xlsx", startRow = 2, colNames = TRUE) 
	normmetab <- subset(metab, metab$Filename %in% metabsamples & metab$Filename != '15316oA_72')
	drop_cols <- c("Genotype","Tissue","Condition","Time.(day)","Replicate","Weight.(mg)")
	normmetab <- normmetab[, !colnames(normmetab) %in% drop_cols]
	rownames(normmetab) <- normmetab$Filename
	normmetab$Filename <- NULL
	normmetab <- t(normmetab)

	# replace metab ids with matched expression ones
	matchids_col <- function(df, oldids, newids) {
		i <- match(oldids, colnames(df), nomatch=0)
		colnames(df)[i] <- as.character(newids)[i]
		return(df)
	}
	
	#place=colnames(df)
	oldids=pairedsamples$Metabolomics
	newids=pairedsamples$RNAseq
	df=normmetab
	normmetabids <- matchids_col(df, oldids, newids)
	save(normmetabids, file=fname)
}else{load(fname)} 


# check metabolite id match between arabidopsis and ccle network
# load ccle network
load('~/Documents/metabpredplants/data/ccle-regulon.rda')
	
fname <- "~/Documents/metabpredplants/LozanoElena2022/data/finalmetabmat.rda"
if(!file.exists(fname)){
	cclemetabolites <- tolower(names(network))
	arabmetabolites <- tolower(rownames(normmetabids))
	#write.table(arabmetabolites, "~/Documents/metabpredplants/LozanoElena2022/data/arabmetabolites_lozano.csv", sep='\t', quote=F, row.names=F)
	commonmetabolites <- intersect(cclemetabolites, arabmetabolites)
	
	normmetabids1 <- subset(normmetabids, tolower(rownames(normmetabids)) %in% commonmetabolites)
	#What is in ccle that is not in arabidopsis?
	# inccle <- setdiff(cclemetabolites, arabmetabolites) 
	#write.table(inccle, "~/Documents/metabpredplants/data/inccle.csv", sep='\t', quote=F, row.names=F)
	#What is in arabidopsis that is not in ccle?
	# inarab <- setdiff(arabmetabolites, cclemetabolites) 
	# write.table(inarab, "~/Documents/metabpredplants/data/inarab.csv", sep='\t', quote=F, row.names=F)
	#ccle					-> arab 
	#"f1p/f6p/g1p/g6p" 	-> "glucose.6.phosphate" / "glucose.1.phosphate"

	# replace metabolite ids recovered by pubchem synonyms
	pubchem <- read.delim("~/Documents/metabpredplants/data/pubchem_syns.csv", sep='\t')
	normmetabpubchem  <- subset(normmetabids, tolower(rownames(normmetabids)) %in% pubchem$arab)
	
	oldids=pubchem$arab
	newids=pubchem$ccle
	df=normmetabpubchem
	
	matchids_row <- function(df, oldids, newids) {
		i <- match(tolower(oldids),tolower(rownames(df)), nomatch=0)
		rownames(df)[i] <- as.character(newids)[i]
		return(df)
		}	

	normmetabids2 <- matchids_row(df, oldids, newids)
	metabmat <- rbind(normmetabids2,normmetabids1)
	metabmatfinal <- metabmat[complete.cases(metabmat), ] # remove NA: Isocitrate
	# first letter down, except for GABA: capital letter

#	newnames=c()	
#	firstUp <- function(df) {
#		for (i in rownames(df)){
#			if (i == toupper(i)){
#				return(i)
#			}else{
#				substr(i, 1, 1) <- toupper(substr(i, 1, 1))
#				newnames = c(newnames,i)
#				rownames(df) <- newnames
#  				return(df)}
# 				}
#  			}
# Sometimes just a 'string' library [laughs]
 			
	rownames(metabmatfinal)[1:5] <- str_to_title(rownames(metabmatfinal)[1:5])
	save(metabmatfinal, file=fname)
}else{load(fname)} 

write.table(metabmatfinal, "metabmatfinal.csv", sep='\t', row.names=T, quote=F)  


## Final matrix 
# create a unique matrix of transcript and metabolite levels
# fname <- "~/Documents/metabpredplants/LozanoElena2022/data/metabexp.rda"

# if(!file.exists(fname)){
# 	metabexp <- rbind(normcounts, metabmatfinal)
# 	save(metabexp, file=fname)
# }else{load(fname)} 


## Metabolite prediction: sample_by_sample=sbs  ----------------------------
pred_lozano <- mra(normcounts, regulon = network)
common <- intersect(rownames(metabmatfinal), rownames(pred_lozano))

# predicted/real correlation
source("~/Documents/metabpredplants/code/pred_real_correlation.R")
thr=0.3
pred_real_correlation <- correlation(metabmatfinal[common,], pred_lozano[common,], thr)
#[1] "pyroglutamic acid" "inositol"          "gaba"             
#[4] "lysine"            "methionine"        "phenylalanine"    
#[7] "putrescine"        "tyrosine"

lozano_corto <- pred_real_correlation[[1]]$cors_spearman #spearman
lozano_cortosign <- pred_real_correlation[[2]]
save(lozano_corto, file = "~/Documents/metabpredplants/LozanoElena2022/results/lozano_corto.rda")

wilcox_lozano <- wilcox.test(lozano_corto, alternative = "greater")
#data:  lozano_corto
#V = 73.5, p-value = 0.7072


#density plot of pred vs real correlations
library(grid)
library(gridExtra)

png("Lozano_densityCC_sbs.png", height=15, width=17, res=600, units="cm")
# plot(density(lozano_corto), lwd=2, xlim=c(-1,1), xlab="Correlation coefficient predicted vs measured", main="", col='royalblue1')
# mtext(paste0("Wilcoxon signed rank test pvalue = ", round(wilcox_lozano$p.value, 3)))
# abline(v=0, lty=2, col="gray", lwd=2)
l <- density(lozano_corto)
p <- ggplot(pred_real_correlation[[1]], aes(x=lozano_corto)) + 
	geom_density(color='black', size=1, fill = "royalblue1", alpha = 0.2) +
	xlim(range(l$x)) +
	theme_bw() +
	xlab("Correlation coefficient predicted vs measured") +
  	ylab("Density") +
  	theme(axis.text=element_text(size=15, color = "black"), axis.title=element_text(size=15))

# Add line
pp <- p + geom_vline(aes(xintercept=0),
	color="grey", linetype="dashed", size=1.5)
#round(wilcox_lozano$p.value, 3)=0.707
grid.arrange(pp, top = textGrob("Wilcoxon signed rank test pvalue = 0.707", rot = 0, vjust = 1, gp=gpar(col="black", fontsize=14)))
dev.off()

# barchart 
png("Lozano_barplotCC_sbs.png", height=15, width=17, res=600, units="cm")
values=pred_real_correlation[[1]]$cors_spearman
labels=pred_real_correlation[[1]]$metab_spearman
  ggplot(pred_real_correlation[[1]], aes(x = reorder(labels, -values), y = values)) +
  geom_bar(stat = "identity", fill='royalblue1') +
  theme_bw() +
  coord_flip() +
  theme(legend.position = "none")  +
  theme(axis.text=element_text(size=15, color = "black"), axis.title=element_text(size=15)) +
  ggtitle("") +
  theme(plot.title = element_text(size=15))+
  xlab("")+
  ylab("Spearman correlation predicted vs measured")   
dev.off()


