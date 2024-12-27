## Metabolite level prediction from gene expression data of Arabidopsis in
## Lozano-Elena et al. 2022 
## by ccle network
## May 24th, 2022                                                                
## *************************************

## Settings ----------------------------
set.seed(1)
rm(list=ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, plyr, readr, tidyverse, fs, biomaRt)	

setwd("~/Documents/metabpredplants/LozanoElena2022")

# load normalised expression data
files <- dir_ls("~/Documents/metabpredplants/LozanoElena2022/data/orthologs", recurse = TRUE, type = "file", glob = "*.csv")  

fname <- "~/Documents/metabpredplants/LozanoElena2022/data/orthologs/orthologs.rda"
		
## DIOPT Ortholog prediction - output parsing ----------------------------

ortholog <- function(files){
	list_diopt <- map(files, read.csv)
	diopt <- (Reduce(function(...) rbind(...), list_diopt))
	diopt <- subset(diopt, diopt$Human.GeneID != "none found") #19434    19
	rankOrder <- c("high", "moderate","low")
	scoreOrder <- c("Yes", "No")
	sorteddiopt<- diopt[order(diopt$Search.Term, -diopt$DIOPT.Score, -diopt$Weighted.Score,  match(diopt$Rank, rankOrder), match(diopt$Best.Score, scoreOrder), match(diopt$Best.Score.Reverse, scoreOrder)), ]
	dioptbestarab <- sorteddiopt[!duplicated(sorteddiopt$Search.Term), ]	
	sorteddioptH <- dioptbestarab[order(dioptbestarab$Human.Symbol, -dioptbestarab$DIOPT.Score, -dioptbestarab$Weighted.Score,  match(dioptbestarab$Rank, rankOrder), match(dioptbestarab$Best.Score, scoreOrder), match(dioptbestarab$Best.Score.Reverse, scoreOrder)), ]
	orthologs <- sorteddioptH[!duplicated(sorteddioptH$Human.Symbol), ]	
	return(orthologs)	
}

#check duplications
#dioptbestH[duplicated(dioptbestH$Search.Term), ]
#dioptbestH[duplicated(dioptbestH$Human.Symbol), ]

# Gene annotation
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

#attributes
#goslim_goa_accession
#goslim_goa_description

#ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
#values = orthologs$Human.Symbol
#annotation <- getBM(attributes=c('hgnc_symbol', 'description' ),filters = 'hgnc_symbol', values = values, mart = ensembl)
#fullOrthologs <- merge(orthologs, annotation, by.x = 'Human.Symbol', by.y = 'hgnc_symbol', all.x=T)

if(!file.exists(fname)){
	orthologs <- ortholog(files)
	save(fullOrthologs, file=fname)
}else{load(fname)} 

write.table(fullOrthologs, "orthologs.csv", sep='\t', row.names=F, quote=F)  
