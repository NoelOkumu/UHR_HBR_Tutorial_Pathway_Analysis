#This is a pathway analysis script for my UHR_HBR data.
#Set the working directory
#Define working dir paths

datadir = "~/RNA_SeqData_Pipeline/RNA_REFS_DIR/RNA_DATA_DIR/Indexing/UHR_HCC_sam_bam/UHC_HCC_bam/HTSeq_UHR_HBR_counts/de/htseq_counts/UHR_HBR_DESeq2_Output"
setwd(datadir)

#Load in the DE results file with only significant genes
DE_genes = read.table("DE_sig_genes_DESeq2.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#Import tools
library(AnnotationDbi)
library(org.Hs.eg.db)
library(GO.db)
library(gage)

# Set up go database
go.hs = go.gsets(species = "human")
go.bp.gs = go.hs$go.sets[go.hs$go.subs$BP]
go.mf.gs = go.hs$go.sets[go.hs$go.subs$MF]
go.cc.gs = go.hs$go.sets[go.hs$go.subs$CC]

# Here I read in an MSigDB gene set that had been selected for the UHR_HBR tutorial data 
c8 = "http://genomedata.org/rnaseq-tutorial/c8.all.v7.2.entrez.gmt"
all_cell_types = readList(c8)

#We now have our DEGs and our Gene sets
head(all_cell_types)

# Read the DE results file with specified column names
DE_genes <- read.table("DE_sig_genes_DESeq2.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")

# Verify the column names
colnames(DE_genes)

DE_genes$entrez = mapIds(org.Hs.eg.db, column = "ENTREZID", keys = DE_genes$ensemblID, keytype = "ENSEMBL", multiVals = "first")

head(DE_genes)

#Remove spike-in
DE_genes_clean = DE_genes[!grepl("ERCC", DE_genes$ensemblID), ]

##Just so we know what we have removed 
ERCC_gene_count = nrow(DE_genes[grepl("ERCC", DE_genes$ensemblID), ])
ERCC_gene_count

###Deal with genes that we do not have an Entrez ID for 
missing_ensembl_key = DE_genes_clean[is.na(DE_genes_clean$entrez), ]
DE_genes_clean = DE_genes_clean[!DE_genes_clean$ensemblID %in% missing_ensembl_key$ensemblID, ]

###Try mapping using a different key
missing_ensembl_key$entrez = mapIds(org.Hs.eg.db, column = "ENTREZID", keys = missing_ensembl_key$Symbol, keytype = "SYMBOL", multiVal = "first")

#Remove remaining genes 
missing_ensembl_key_update = missing_ensembl_key[!is.na(missing_ensembl_key$entrez),]

#Create a Final Gene list of all genes where we were able to find an Entrez ID (using two approaches)
DE_genes_clean = rbind(DE_genes_clean, missing_ensembl_key_update)
head(DE_genes_clean)

# grab the log fold changes for everything
De_gene.fc = DE_genes_clean$log2FoldChange

# set the name for each row to be the Entrez Gene ID
names(De_gene.fc) = DE_genes_clean$entrez

#Run GAGE
#go 
fc.go.bp.p = gage(De_gene.fc, gsets = go.bp.gs)
fc.go.mf.p = gage(De_gene.fc, gsets = go.mf.gs)
fc.go.cc.p = gage(De_gene.fc, gsets = go.cc.gs)

#msigdb
fc.c8.p = gage(De_gene.fc, gsets = all_cell_types)

###Convert to dataframes 
#Results for testing for GO terms which are up-regulated
fc.go.bp.p.up = as.data.frame(fc.go.bp.p$greater)
fc.go.mf.p.up = as.data.frame(fc.go.mf.p$greater)
fc.go.cc.p.up = as.data.frame(fc.go.cc.p$greater)

#Results for testing for GO terms which are down-regulated
fc.go.bp.p.down = as.data.frame(fc.go.bp.p$less)
fc.go.mf.p.down = as.data.frame(fc.go.mf.p$less)
fc.go.cc.p.down = as.data.frame(fc.go.cc.p$less)

#Results for testing for MSigDB C8 gene sets which are up-regulated
fc.c8.p.up = as.data.frame(fc.c8.p$greater)
head(fc.c8.p.up)

#Results for testing for MSigDB C8 gene sets which are down-regulated
fc.c8.p.down = as.data.frame(fc.c8.p$less)
head(fc.c8.p.down)

#View the top 20 significantly up- or down-regulated GO terms from the Cellular Component Ontology
head(fc.go.cc.p.up[order(fc.go.cc.p.up$p.val),], n = 20)
head(fc.go.cc.p.down[order(fc.go.cc.p.down$p.val),], n = 20)

#What if we want to know which specific genes from our DE gene result were found in a specific significant pathway?
#For example, one significant pathway from fc.go.cc.p.down was "GO:0098794 postsynapse" with set.size = 11 genes.
#Let's extract the postsynapse DE gene results
postsynapse = DE_genes_clean[which(DE_genes_clean$entrez %in% go.cc.gs$`GO:0098794 postsynapse`),]

#How many total postsynapse genes are there in GO? How many total DE genes? How many overlap?
length(go.cc.gs$`GO:0098794 postsynapse`)
length(DE_genes_clean$entrez)
length(postsynapse)

#Are the postsynapse DE genes consistently down-regulated? Let's print out a subset of columns from the DE result for postsynapse genes
postsynapse[,c("Symbol", "entrez", "log2FoldChange", "padj", "UHR_Rep1", "UHR_Rep2", "UHR_Rep3", "HBR_Rep1", "HBR_Rep2", "HBR_Rep3")]

write.table(fc.go.cc.p.up, "fc.go.cc.p.up.tsv", quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
write.table(fc.go.cc.p.down, "fc.go.cc.p.down.tsv", quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)

list.files()

quit(save = "no")

getwd()

list.files()




