library(DESeq2)
library(stringr)
library(sqldf)
library(reshape2)
library(ggplot2)
library(gridExtra)

# 
# plot_volcano <- function(res, pcutoff=200,cex=0.5){
#   #res <- merge(res, gene_anno)
#   logp <- -log10(res$pvalue)
#   logp <- pmin(logp,pcutoff)
#   plot(res$log2FoldChange,logp,cex=0)
#   text(res$log2FoldChange,logp, labels = res$symbol,cex=cex)
# }


# 
# ##################################
# ## Rename colums before storing deseq comparison
# store_de_list <- function(f,dat){
#   dat <- dat[,c("geneid","log2FoldChange","lfcSE","pvalue","padj","symbol")]
#   dat <- dat[order(dat$pvalue),]
#   write.csv(dat, file = f, row.names = FALSE)
# }

###################
################### Load data
###################

counts <- read.table("raw_rnaseq/counts.csv", header = TRUE, stringsAsFactors = FALSE)
rownames(counts) <- counts$Geneid
counts <- counts[,-(1:6)]

cond <- read.csv("raw_rnaseq/cond.csv", stringsAsFactors = FALSE, row.names = 1)
cond$treatment <- factor(cond$treatment)
cond$sample <- rownames(cond)

## Normalized reads
dds <- DESeqDataSetFromMatrix(countData = counts, colData = cond, design = ~1)
dds <- DESeq(dds) #overkill
cnt_norm <- counts(dds, normalize=TRUE)

### Prepare gene annotation. was made by:
#	cut -f 1,2,3 fbgn_annotation_ID_fb_2020_03.tsv > gene_annot.csv
gene_anno <- read.table("gene_annot.csv",stringsAsFactors = FALSE)
colnames(gene_anno) <- c("symbol","organism","geneid")
gene_anno$recep_or <- str_sub(gene_anno$symbol,1,2)=="Or" & 
  !(str_sub(gene_anno$symbol,1,3)=="Orc") & 
  !(str_sub(gene_anno$symbol,1,3)=="Ork") & 
  !(str_sub(gene_anno$symbol,1,3)=="Ora")
gene_anno[gene_anno$symbol=="Orco",]$recep_or <- TRUE
gene_anno$recep_or_orig <- gene_anno$recep_or

gene_anno$recep_ir <- str_sub(gene_anno$symbol,1,2)=="Ir" & 
  !(str_sub(gene_anno$symbol,1,3)=="Ire") & 
  !(str_sub(gene_anno$symbol,1,3)=="Irb") & 
  !(str_sub(gene_anno$symbol,1,3)=="Irk") & 
  !(str_sub(gene_anno$symbol,1,3)=="Irc") & 
  !(str_sub(gene_anno$symbol,1,3)=="Irp") & 
  !(str_sub(gene_anno$symbol,1,3)=="Iri") & 
  !(str_sub(gene_anno$symbol,1,3)=="Iru")

gene_anno$recep_gr <- str_sub(gene_anno$symbol,1,2)=="Gr" & 
  !(str_sub(gene_anno$symbol,1,3)=="Gra") & 
  !(str_sub(gene_anno$symbol,1,3)=="Grd") & 
  !(str_sub(gene_anno$symbol,1,3)=="Gri") & 
  !(str_sub(gene_anno$symbol,1,3)=="Grx") 

gene_anno$other_int_che <- str_sub(gene_anno$symbol,1,3)=="Che"
gene_anno$other_int_obp <- str_sub(gene_anno$symbol,1,3)=="Obp"

#### Make a class variable from all previous groups
gene_anno$gclass <- ""
gene_anno$gclass[gene_anno$recep_ir] <- "Ir"
gene_anno$gclass[gene_anno$recep_or] <- "Or"
gene_anno$gclass[gene_anno$recep_gr] <- "Gr"
gene_anno$gclass[gene_anno$symbol=="Orco"] <- "Orco"
gene_anno$gclass[gene_anno$other_int_che] <- "Che"
gene_anno$gclass[gene_anno$other_int_obp] <- "Obp"

#geneid_oe <- gene_anno$geneid[gene_anno=="Or65a"]  #which?

#sym_exp_receptors_or <- gene_anno$symbol[gene_anno$gclass=="Or"]
sym_exp_receptors_or <- read.table("list_or.csv",stringsAsFactors = FALSE)[,1]
gene_anno$recep_or <- gene_anno$symbol %in% sym_exp_receptors_or
sym_exp_receptors_ir <- gene_anno$symbol[gene_anno$gclass=="Ir"]
sym_exp_receptors_gr <- gene_anno$symbol[gene_anno$gclass=="Gr"]

genesym_receptors <- c(sym_exp_receptors_or, sym_exp_receptors_ir, sym_exp_receptors_gr)
geneid_receptors <- gene_anno$geneid[gene_anno$symbol %in% genesym_receptors]






########################## should refactor this code!!!


#### Other genes
intgene <- read.csv("genelist.csv", stringsAsFactors = FALSE)
intgene$isin <- intgene$symbol %in% gene_anno$symbol

gene_anno$is_housekeeping <- gene_anno$symbol %in% intgene$symbol

gene_anno$is_synapse_osn <- gene_anno$symbol %in% intgene$symbol[intgene$category %in% c("newsynapse","olfactory sensory neuron markers")]
gene_anno$is_cilia <- gene_anno$symbol %in% intgene$symbol[intgene$category %in% c("cilia")]

gene_anno$is_trichoids <- gene_anno$symbol %in% intgene$symbol[intgene$category %in% c("trichoids")]
gene_anno$is_basiconic <- gene_anno$symbol %in% intgene$symbol[intgene$category %in% c("basiconic")]



