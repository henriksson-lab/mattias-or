### orco bort i orco. markörgener med elav. 1-4 markörer
### or42b,  not or65a!


library(DESeq2)
library(stringr)
library(sqldf)
library(reshape2)
library(ggplot2)
library(gridExtra)


plot_volcano <- function(res, pcutoff=200,cex=0.5){
  #res <- merge(res, gene_anno)
  logp <- -log10(res$pvalue)
  logp <- pmin(logp,pcutoff)
  plot(res$log2FoldChange,logp,cex=0)
  text(res$log2FoldChange,logp, labels = res$symbol,cex=cex)
}

store_de_list <- function(f,dat){
  dat <- dat[,c("geneid","log2FoldChange","lfcSE","pvalue","padj","symbol")]
  dat <- dat[order(dat$pvalue),]
  write.csv(dat, file = f, row.names = FALSE)
}

###################
################### Load data
###################

counts <- read.table("counts.csv", header = TRUE, stringsAsFactors = FALSE)
rownames(counts) <- counts$Geneid
counts <- counts[,-(1:6)]

cond <- read.csv("cond.csv", stringsAsFactors = FALSE, row.names = 1)
cond$treatment <- factor(cond$treatment)
cond$sample <- rownames(cond)

## Normalized reads
dds <- DESeqDataSetFromMatrix(countData = counts, colData = cond, design = ~1)
dds <- DESeq(dds) #overkill
cnt_norm <- counts(dds, normalize=TRUE)

### Prepare gene annotation
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
#gene_anno$symbol[gene_anno$other_int]
#pick above 1.3

gene_anno$gclass <- ""
gene_anno$gclass[gene_anno$recep_ir] <- "Ir"
gene_anno$gclass[gene_anno$recep_or] <- "Or"
gene_anno$gclass[gene_anno$recep_gr] <- "Gr"
gene_anno$gclass[gene_anno$symbol=="Orco"] <- "Orco"
#gene_anno$gclass[gene_anno$symbol=="Or65a"] <- "Or65a"
#gene_anno$gclass[gene_anno$symbol=="Or47b"] <- "Or47b"
gene_anno$gclass[gene_anno$other_int_che] <- "Che"
gene_anno$gclass[gene_anno$other_int_obp] <- "Obp"
table(gene_anno$gclass)

geneid_oe <- gene_anno$geneid[gene_anno=="Or65a"]  #which?

#sym_exp_receptors_or <- gene_anno$symbol[gene_anno$gclass=="Or"]
sym_exp_receptors_or <- read.table("list_or.csv",stringsAsFactors = FALSE)[,1]
gene_anno$recep_or <- gene_anno$symbol %in% sym_exp_receptors_or
sym_exp_receptors_ir <- gene_anno$symbol[gene_anno$gclass=="Ir"]
sym_exp_receptors_gr <- gene_anno$symbol[gene_anno$gclass=="Gr"]

genesym_receptors <- c(sym_exp_receptors_or, sym_exp_receptors_ir, sym_exp_receptors_gr)
geneid_receptors <- gene_anno$geneid[gene_anno$symbol %in% genesym_receptors]

#### Other genes
intgene <- read.csv("genelist.csv", stringsAsFactors = FALSE)
intgene$isin <- intgene$symbol %in% gene_anno$symbol
#intgene[!intgene$isin,]

gene_anno$is_housekeeping <- gene_anno$symbol %in% intgene$symbol#[intgene$category %in% c("housekeepingExp","synapse markers")]

gene_anno$is_synapse_osn <- gene_anno$symbol %in% intgene$symbol[intgene$category %in% c("newsynapse","olfactory sensory neuron markers")]
gene_anno$is_cilia <- gene_anno$symbol %in% intgene$symbol[intgene$category %in% c("cilia")]

gene_anno$is_trichoids <- gene_anno$symbol %in% intgene$symbol[intgene$category %in% c("trichoids")]
gene_anno$is_basiconic <- gene_anno$symbol %in% intgene$symbol[intgene$category %in% c("basiconic")]

###################
################### Overall DE: OE vs WT (not used)
###################

dds <- DESeqDataSetFromMatrix(countData = counts, colData = cond, design = ~ day+treatment)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast = c("treatment","oe","wt"))
res <- as.data.frame(res)
res$geneid <- rownames(res)

res <- merge(res, gene_anno)
res <- res[order(res$padj),]

res_oe_orig <- res
res_oe <- data.frame(
  symbol=res$symbol,
  oe_fc=res$log2FoldChange,
  oe_padj=res$padj)


###################
################### overall DE: Compare octo-mutant vs WT
###################

keepdata <- cond$day=="14"
dds <- DESeqDataSetFromMatrix(countData = counts[,keepdata], colData = cond[keepdata,], design = ~ treatment)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast = c("treatment","orco","wt"))
res <- as.data.frame(res)
res$geneid <- rownames(res)

res <- merge(res, gene_anno)
res <- res[order(res$padj),]

res_orco_orig <- res
res_orco <- data.frame(
  symbol=res$symbol,
  orco_fc=res$log2FoldChange,
  orco_padj=res$padj)




###################
################### Histogram of exp levels
###################

keepdata <- cond$day==4 & cond$treatment=="wt"
dds <- DESeqDataSetFromMatrix(countData = counts[,keepdata], colData = cond[keepdata,], design = ~ 1)
dds <- DESeq(dds)
res <- results(dds)#, contrast = c("treatment","oe","wt"))
res <- as.data.frame(res)
res$geneid <- rownames(res)
res <- merge(res, gene_anno)
res <- res[order(res$baseMean),]
res$symbol <- factor(res$symbol, levels = res$symbol)

df <- data.frame(x=log10(1+res$baseMean), symbol=res$symbol)
df$is_or <- df$symbol %in% sym_exp_receptors_or
df$is_gr <- df$symbol %in% sym_exp_receptors_gr
df$is_ir <- df$symbol %in% sym_exp_receptors_ir

### one method
p1 <- ggplot(df, aes(x=x)) +  
  geom_histogram(binwidth=0.1) + xlim(0,5) + xlab("Overall gene levels")
p2 <- ggplot(df[df$is_or,], aes(x=x)) +  
  geom_histogram(binwidth=0.1) + xlim(0,5) + xlab("Or*")
p3 <- ggplot(df[df$is_ir,], aes(x=x)) +  
  geom_histogram(binwidth=0.1) + xlim(0,5) + xlab("Ir*")
p4 <- ggplot(df[df$is_gr,], aes(x=x)) +  
  geom_histogram(binwidth=0.1) + xlim(0,5) + xlab("Gr*")
gc <- grid.arrange(p1, p2, p3, p4, nrow = 2)


## second method
p1 <- ggplot(df, aes(x=x)) +  
  geom_histogram(binwidth=0.1) + xlim(0,5) + xlab("Overall gene levels")
p2<-ggplot(data=df[df$is_or,], aes(x=symbol, y=x)) +
  geom_bar(stat="identity") + 
  coord_flip() + ylim(0,5)
p3<-ggplot(data=df[df$is_ir,], aes(x=symbol, y=x)) +
  geom_bar(stat="identity") + 
  coord_flip() + ylim(0,5)
p4<-ggplot(data=df[df$is_gr,], aes(x=symbol, y=x)) +
  geom_bar(stat="identity") + 
  coord_flip() + ylim(0,5)
gc <- grid.arrange(p1, p2, p3, p4, nrow = 2)
gc
ggsave("out/genelevels.pdf",gc, width = 10, height = 13)

 


###################
################### Plots of wildtype
###################

#CPM normalize
rec_count <- cnt_norm[,cond$treatment=="wt"]
for(i in 1:ncol(rec_count)){
  rec_count[,i] <- rec_count[,i]/(sum(rec_count[,i])/1e6)
}

rec_count <- log10(1+rec_count)
#rec_count <- log10(1+cnt_norm[rownames(cnt_norm) %in% geneid_receptors,cond$treatment=="wt"])
rec_count <- melt(rec_count)
colnames(rec_count) <- c("geneid","sample","count")
rec_count <- merge(rec_count, cond)[,c("geneid","count","day","treatment")]

rec_count_proc <- sqldf("select avg(count) as count, stdev(count) as var, geneid, day from rec_count group by geneid, day")
rec_count_proc <- merge(rec_count_proc, gene_anno)
rec_count_proc <- rec_count_proc[order(rec_count_proc$count, decreasing = TRUE),]

rec_count_proc <- merge(rec_count_proc,sqldf("select count as count1, geneid from rec_count_proc where day=1"))
rec_count_proc <- merge(rec_count_proc,sqldf("select count as count4, geneid from rec_count_proc where day=4"))
rec_count_proc <- merge(rec_count_proc,sqldf("select count as count14, geneid from rec_count_proc where day=14"))
rec_count_proc$isinc <- rec_count_proc$count14 > rec_count_proc$count4
rec_count_proc$isinc14_4 <- rec_count_proc$count14 > rec_count_proc$count4
rec_count_proc$isinc4_1 <- rec_count_proc$count4 > rec_count_proc$count1



#### Decide which genes are "expressed"
rec_count_proc$isexp <- TRUE
cutoff_gr <- rec_count_proc[rec_count_proc$day == 4 & rec_count_proc$symbol=="Gr28a",]$count4
cutoff_ir <- rec_count_proc[rec_count_proc$day == 4 & rec_count_proc$symbol=="Ir60a",]$count4
cutoff_or <- rec_count_proc[rec_count_proc$day == 4 & rec_count_proc$symbol=="Or65c",]$count4   #prevexp <- read.table("expGenes.txt", stringsAsFactors = FALSE)[,1]   
rec_count_proc$isexp[rec_count_proc$recep_gr] <- rec_count_proc$count4[rec_count_proc$recep_gr] >= cutoff_gr
rec_count_proc$isexp[rec_count_proc$recep_ir] <- rec_count_proc$count4[rec_count_proc$recep_ir] >= cutoff_ir
#rec_count_proc$isexp[rec_count_proc$recep_or] <- rec_count_proc$count4[rec_count_proc$recep_or] >= cutoff_or

#### Transfer back to symbol table
gene_anno$isexp <- gene_anno$symbol %in% rec_count_proc$symbol[rec_count_proc$isexp]

rec_count_proc$isexp_basal <- rec_count_proc$count4 > 1.3


################## PLOTXY
################## PLOTXY OR GR IR
################## PLOTXY

rec_count_proc$colcat <- "gray"
rec_count_proc$colcat[rec_count_proc$recep_ir & rec_count_proc$isexp] <- "red"
rec_count_proc$colcat[rec_count_proc$recep_gr & rec_count_proc$isexp] <- "blue"
rec_count_proc$colcat[rec_count_proc$recep_or & rec_count_proc$isexp] <- "green"
# rec_count_proc$colcat[rec_count_proc$recep_ir] <- "red"
# rec_count_proc$colcat[rec_count_proc$recep_gr] <- "blue"
rec_count_proc$colcat[rec_count_proc$recep_or_orig] <- "green"

plotxy("count1","count4","Log10(WT day 1 CPM)","Log10(WT day 4 CPM)","newerout/xy_wt_1vs4_orgrir.pdf")
plotxy("count4","count14","Log10(WT day 4 CPM)","Log10(WT day 14 CPM)","newerout/xy_wt_4vs14_orgrir.pdf")

############### Obp

rec_count_proc$colcat <- "gray"
rec_count_proc$colcat[rec_count_proc$other_int_obp & rec_count_proc$isexp_basal] <- "red"
#rec_count_proc$colcat[rec_count_proc$other_int_che & rec_count_proc$isexp_basal] <- "blue"

plotxy("count1","count4","Log10(WT day 1 CPM)","Log10(WT day 4 CPM)","newerout/xy_wt_1vs4_obp.pdf")
plotxy("count4","count14","Log10(WT day 4 CPM)","Log10(WT day 14 CPM)","newerout/xy_wt_4vs14_obp.pdf")

############### synapse and OSN markers

rec_count_proc$colcat <- "gray"
rec_count_proc$colcat[rec_count_proc$is_synapse_osn] <- "red"  #& rec_count_proc$isexp_basal
rec_count_proc$colcat[rec_count_proc$symbol=="elav"] <- "blue"

plotxy("count1","count4","Log10(WT day 1 CPM)","Log10(WT day 4 CPM)","newerout/xy_wt_1vs4_synapse.pdf")
plotxy("count4","count14","Log10(WT day 4 CPM)","Log10(WT day 14 CPM)","newerout/xy_wt_4vs14_synapse.pdf")

############### ciliome

rec_count_proc$colcat <- "gray"
rec_count_proc$colcat[rec_count_proc$is_cilia] <- "red"

plotxy("count1","count4","Log10(WT day 1 CPM)","Log10(WT day 4 CPM)","newerout/xy_wt_1vs4_cilia.pdf")
plotxy("count4","count14","Log10(WT day 4 CPM)","Log10(WT day 14 CPM)","newerout/xy_wt_4vs14_cilia.pdf")


############### d1 vs d4

rec_count_proc$colcat <- "gray"
rec_count_proc$colcat[rec_count_proc$symbol %in% intgene$symbol[intgene$category=="listAntenna"]] <- "red"
plotxy("count1","count4","Log10(WT day 1 CPM)","Log10(WT day 4 CPM)","out.20201015/xy_wt_1vs4_antenna.pdf")

rec_count_proc$colcat <- "gray"
rec_count_proc$colcat[rec_count_proc$symbol %in% intgene$symbol[intgene$category=="listLAV"]] <- "red"
plotxy("count1","count4","Log10(WT day 1 CPM)","Log10(WT day 4 CPM)","out.20201015/xy_wt_1vs4_LAV.pdf")

rec_count_proc$colcat <- "gray"
rec_count_proc$colcat[rec_count_proc$symbol %in% intgene$symbol[intgene$category=="listLE"]] <- "red"
plotxy("count1","count4","Log10(WT day 1 CPM)","Log10(WT day 4 CPM)","out.20201015/xy_wt_1vs4_LE.pdf")



################## Time course plots, not separated by level
p1 <- ggplot(rec_count_proc[rec_count_proc$gclass=="Ir",], aes(x = day, y = count, colour = isexp, group=geneid)) + geom_line() 
p2 <- ggplot(rec_count_proc[rec_count_proc$gclass=="Or",], aes(x = day, y = count, colour = prevexp, group=geneid)) + geom_line() 
p3 <- ggplot(rec_count_proc[rec_count_proc$gclass=="Gr",], aes(x = day, y = count, colour = isexp, group=geneid)) + geom_line() 
gc <- grid.arrange(p1, p2, p3, nrow = 1)

################### Time course plot, going up vs going down
p1 <- ggplot(rec_count_proc[rec_count_proc$isinc14_4 & rec_count_proc$isinc4_1,], aes(x = day, y = count, colour = gclass, group=geneid)) + 
  geom_line() + ylim(0,5)
p2 <- ggplot(rec_count_proc[!rec_count_proc$isinc14_4 & rec_count_proc$isinc4_1,], aes(x = day, y = count, colour = gclass, group=geneid)) + 
  geom_line() + ylim(0,5)
p3 <- ggplot(rec_count_proc[rec_count_proc$isinc14_4 & !rec_count_proc$isinc4_1,], aes(x = day, y = count, colour = gclass, group=geneid)) + 
  geom_line()  + ylim(0,5)
p4 <- ggplot(rec_count_proc[!rec_count_proc$isinc14_4 & !rec_count_proc$isinc4_1,], aes(x = day, y = count, colour = gclass, group=geneid)) + 
  geom_line() + ylim(0,5)
gc <- grid.arrange(p1, p2, p3, p4, nrow = 2)
#ggsave("out/tc_wt_4.pdf",gc)



###################
################### OE time course (not used)
###################

rec_count <- log10(1+cnt_norm[rownames(cnt_norm) %in% geneid_receptors,cond$treatment=="oe"])
rec_count <- melt(rec_count)
colnames(rec_count) <- c("geneid","sample","count")
rec_count <- merge(rec_count, cond)[,c("geneid","count","day","treatment")]

rec_count_proc <- sqldf("select avg(count) as count, stdev(count) as var, geneid, day from rec_count group by geneid, day")
rec_count_proc <- merge(rec_count_proc, gene_anno)
rec_count_proc <- rec_count_proc[order(rec_count_proc$count, decreasing = TRUE),]

rec_count_proc <- merge(rec_count_proc,sqldf("select count as count1, geneid from rec_count_proc where day=1"))
rec_count_proc <- merge(rec_count_proc,sqldf("select count as count4, geneid from rec_count_proc where day=4"))
rec_count_proc$isinc <- rec_count_proc$count4 > rec_count_proc$count1

p <- ggplot(rec_count_proc, aes(x = day, y = count, colour = gclass)) + 
  geom_line() +
  theme(legend.position = "none") 
p

p1 <- ggplot(rec_count_proc[rec_count_proc$gclass=="Ir",], aes(x = day, y = count, colour = gclass, group=geneid)) + geom_line() 
p2 <- ggplot(rec_count_proc[rec_count_proc$gclass=="Or",], aes(x = day, y = count, colour = gclass, group=geneid)) + geom_line() 
p3 <- ggplot(rec_count_proc[rec_count_proc$gclass=="Gr",], aes(x = day, y = count, colour = gclass, group=geneid)) + geom_line() 
gc <- grid.arrange(p1, p2, p3, nrow = 1)

p1 <- ggplot(rec_count_proc[rec_count_proc$isinc,], aes(x = day, y = count, colour = gclass, group=geneid)) + 
  geom_line()  + ylim(0,5)
p2 <- ggplot(rec_count_proc[!rec_count_proc$isinc,], aes(x = day, y = count, colour = gclass, group=geneid)) + 
  geom_line()  + ylim(0,5)
gc <- grid.arrange(p1, p2, nrow = 1)
#ggsave("out/tc_oe.pdf",gc)


###################
################### 2d plot: basal control level vs FC orco
###################

dat <- merge(merge(res_orco, recep_baseline_14), gene_anno)
dat

pdf("out/orco vs wt14.pdf")
plot(dat$orco_fc, dat$baseline,cex=0.1, col="gray", pch=19, xlab="orco FC",ylab="Normalized counts wt day 14")
text(dat$orco_fc, dat$baseline, labels = dat$symbol,cex=0.8, col="blue")
dev.off()

dat$logp <- -log10(dat$orco_padj)
dat$logp[dat$orco_fc<0] <- -dat$logp[dat$orco_fc<0]
plot(dat$logp, dat$baseline,cex=0.1, col="gray", pch=19)
text(dat$logp, dat$baseline, labels = dat$symbol,cex=0.8, col="blue")



###################
################### volcanoes
###################

plot_volcano(res_orco_orig, pcutoff = 30)
pdf("out/volcano_orco.pdf")
plot_volcano(res_orco_orig, pcutoff = 30)
dev.off()

plot_volcano(res_oe_orig, pcutoff = 50, cex=0.8)
pdf("out/volcano_oe.pdf")
plot_volcano(res_oe_orig, pcutoff = 50, cex=0.8)
dev.off()

plot_volcano(res_orco_orig[res_orco_orig$symbol %in% sym_receptors,], pcutoff = 100,cex=0.8)
pdf("out/volcano_orco_recep.pdf")
plot_volcano(res_orco_orig[res_orco_orig$symbol %in% sym_receptors,], pcutoff = 100,cex=0.8)
dev.off()




###################
################### Barplots of fold-change
###################

orderggByExp <- function(res){
  res$symbol <- factor(res$symbol, levels = c(setdiff(res$symbol,"Orco"), "Orco"))
  res  
}


######################################################### WT day 1 vs 4
cond$dday <- factor(sprintf("d%d",cond$day))
keepdata <- cond$treatment=="wt" & cond$day %in% c(1,4)
dds <- DESeqDataSetFromMatrix(countData = counts[,keepdata], colData = cond[keepdata,], design = ~ dday)
#keep <- rowSums(counts(dds)) >= 10
#dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast = c("dday","d4","d1"))
res <- as.data.frame(res)
res$geneid <- rownames(res)
res <- merge(res, gene_anno)
#res <- res[res$symbol %in% genesym_receptors,]
res_wt14 <- res

### want the d4 here as well
w <- cnt_norm[,cond$treatment=="wt" & cond$day==4]
for(i in 1:ncol(w)){
  w[,i] <- w[,i]/(sum(w[,i])/1e6)
}
w <- log10(1+w)
w <- melt(w)
colnames(w) <- c("geneid","sample","count")
w <- merge(w, cond)#[,c("geneid","count","day","treatment")]
w <- sqldf("select avg(count) as d4count, geneid from w group by geneid")
res_wt14 <- merge(res_wt14,w)

res_wt14 <- res_wt14[order(res_wt14$d4count, decreasing = FALSE),]
res_wt14 <- orderggByExp(res_wt14)
store_de_list("de/wt_d1_vs_d4.csv",res_wt14)


######################################################### WT day 4 vs 14
cond$dday <- factor(sprintf("d%d",cond$day))
keepdata <- cond$treatment=="wt" & cond$day %in% c(4,14)
dds <- DESeqDataSetFromMatrix(countData = counts[,keepdata], colData = cond[keepdata,], design = ~ dday)
#keep <- rowSums(counts(dds)) >= 10
#dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast = c("dday","d14","d4"))
res <- as.data.frame(res)
res$geneid <- rownames(res)
res <- merge(res, gene_anno)
#res <- res[res$symbol %in% genesym_receptors,]
res <- res[order(res$baseMean, decreasing = FALSE),]
res <- orderggByExp(res)
res_wt414 <- res
store_de_list("de/wt_d14_vs_d4.csv",res_wt414)


######################################################### OE day 4 vs 1  (Or65a)
cond$dday <- factor(sprintf("d%d",cond$day))
keepdata <- cond$treatment=="oe" & cond$day %in% c(4,1)
dds <- DESeqDataSetFromMatrix(countData = counts[,keepdata], colData = cond[keepdata,], design = ~ dday)
#keep <- rowSums(counts(dds)) >= 10
#dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast = c("dday","d4","d1"))
res <- as.data.frame(res)
res$geneid <- rownames(res)
res <- merge(res, gene_anno)
#res <- res[res$symbol %in% genesym_receptors,]
res <- res[order(res$baseMean, decreasing = FALSE),]
res <- orderggByExp(res)
res_oe41 <- res
store_de_list("de/oe_d4_vs_d1.csv",res_oe41)


######################################################### OE day 1 vs WT day 1  (Or65a)
cond$dday <- factor(sprintf("d%d",cond$day))
keepdata <- cond$day %in% 1
dds <- DESeqDataSetFromMatrix(countData = counts[,keepdata], colData = cond[keepdata,], design = ~ treatment)
#keep <- rowSums(counts(dds)) >= 10
#dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast = c("treatment","oe","wt"))
res <- as.data.frame(res)
res$geneid <- rownames(res)
res <- merge(res, gene_anno)
#res <- res[res$symbol %in% genesym_receptors,]
res <- res[order(res$baseMean, decreasing = FALSE),]
res <- orderggByExp(res)
res_oe_wt_d1 <- res
store_de_list("de/oe_d1_vs_wt_d1.csv",res_oe_wt_d1)



######################################################### OE day 4 vs WT day 4  (Or65a)
cond$dday <- factor(sprintf("d%d",cond$day))
keepdata <- cond$day %in% 4
dds <- DESeqDataSetFromMatrix(countData = counts[,keepdata], colData = cond[keepdata,], design = ~ treatment)
#keep <- rowSums(counts(dds)) >= 10
#dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast = c("treatment","oe","wt"))
res <- as.data.frame(res)
res$geneid <- rownames(res)
res <- merge(res, gene_anno)
#res <- res[res$symbol %in% genesym_receptors,]
res <- res[order(res$baseMean, decreasing = FALSE),]
res <- orderggByExp(res)
res_oe_wt_d4 <- res
store_de_list("de/oe_d4_vs_wt_d4.csv",res_oe_wt_d4)



######################################################### Orco vs wt
cond$dday <- factor(sprintf("d%d",cond$day))
keepdata <- cond$day %in% 14  #  cond$treatment=="orco" # & cond$day %in% c(4,1)
dds <- DESeqDataSetFromMatrix(countData = counts[,keepdata], colData = cond[keepdata,], design = ~ treatment)
#keep <- rowSums(counts(dds)) >= 10
#dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast = c("treatment","orco","wt"))
res <- as.data.frame(res)
res$geneid <- rownames(res)
res <- merge(res, gene_anno)
#res <- res[res$symbol %in% genesym_receptors,]
res <- res[order(res$baseMean, decreasing = FALSE),]
res <- orderggByExp(res)
res_orco_wt <- res
store_de_list("de/orco_d14_vs_wt_d14.csv",res_orco_wt)




######################################################### The actual plots

### For WT OR
theylim_or <- ylim(-3,3)
p1<-ggplot(data=res_wt14[res_wt14$recep_or,], aes(x=symbol, y=log2FoldChange)) +
  geom_bar(stat="identity") + 
  theylim_or +
  xlab("WT d4 vs d1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2<-ggplot(data=res_wt414[res_wt414$isexp & res_wt414$recep_or,], aes(x=symbol, y=log2FoldChange)) +
  geom_bar(stat="identity") + 
  theylim_or +
  xlab("WT d14 vs d4") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p3<-ggplot(data=res_wt14[res_wt14$recep_or,], aes(x=symbol, y=d4count)) +
  geom_bar(stat="identity") + 
  ylab("log10[WT day 4 CPM]") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p4<-ggplot(data=res_wt414[res_wt414$is_housekeeping,], aes(x=symbol, y=log2FoldChange)) +
#   geom_bar(stat="identity") + 
#   theylim_or +
#   xlab("WT d14 vs d4") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p4<-ggplot(data=res_wt14[res_wt14$is_housekeeping,], aes(x=symbol, y=log2FoldChange)) +
#   geom_bar(stat="identity") +
#   theylim_or +
#   xlab("WT d4 vs d1") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
gc <- grid.arrange(p1, p2,p3, nrow = 2)
ggsave("newout/wt_or.pdf",gc, width = 15, height = 12)

p1<-ggplot(data=res_wt14[res_wt14$symbol %in% intgene$symbol[intgene$category %in% c("olfactory sensory neuron markers")],], aes(x=symbol, y=log2FoldChange)) +
  geom_bar(stat="identity") + 
  theylim_or +
  xlab("WT d4 vs d1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2<-ggplot(data=res_wt14[res_wt14$symbol %in% intgene$symbol[intgene$category %in% c("synapse markers")],], aes(x=symbol, y=log2FoldChange)) +
  geom_bar(stat="identity") + 
  theylim_or +
  xlab("WT d4 vs d1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p3<-ggplot(data=res_wt414[res_wt414$symbol %in% intgene$symbol[intgene$category %in% c("olfactory sensory neuron markers")],], aes(x=symbol, y=log2FoldChange)) +
  geom_bar(stat="identity") + 
  theylim_or +
  xlab("WT d14 vs d4") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p4<-ggplot(data=res_wt414[res_wt414$symbol %in% intgene$symbol[intgene$category %in% c("synapse markers")],], aes(x=symbol, y=log2FoldChange)) +
  geom_bar(stat="identity") + 
  theylim_or +
  xlab("WT d14 vs d4") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
gc <- grid.arrange(p1, p2,p3, p4, nrow = 2)
ggsave("newout/wt_or_ctrl.pdf",gc, width = 15, height = 12)


### For WT IR/GR
theylim_ir <- ylim(-2.5,2.5)
theylim_gr <- ylim(-2.5,2.5)
p1<-ggplot(data=res_wt14[res_wt14$isexp & res_wt14$recep_ir,], aes(x=symbol, y=log2FoldChange)) +
  geom_bar(stat="identity") + 
  theylim_ir +
  xlab("WT d4 vs d1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2<-ggplot(data=res_wt414[res_wt414$isexp & res_wt414$recep_ir,], aes(x=symbol, y=log2FoldChange)) +
  geom_bar(stat="identity") + 
  theylim_ir +
  xlab("WT d14 vs d4") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
gc <- grid.arrange(p1, p2, nrow = 1)
ggsave("newout/wt_ir.pdf",gc, width = 15, height = 8)

### For WT day 1 vs 4, gr
p3<-ggplot(data=res_wt14[res_wt14$isexp & res_wt14$recep_gr,], aes(x=symbol, y=log2FoldChange)) +
  geom_bar(stat="identity") + 
  theylim_gr +
  xlab("WT d4 vs d1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p4<-ggplot(data=res_wt414[res_wt414$isexp & res_wt414$recep_gr,], aes(x=symbol, y=log2FoldChange)) +
  geom_bar(stat="identity") + 
  theylim_gr +
  xlab("WT d14 vs d4") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
gc <- grid.arrange(p1, p2, p3, p4, nrow = 2)
ggsave("newout/wt_gr.pdf",gc, width = 15, height = 8)


### For OE Or65a d4 vs d1 -- OR
p1<-ggplot(data=res_oe41[res_oe41$isexp & res_oe41$recep_or,], aes(x=symbol, y=log2FoldChange)) +
  geom_bar(stat="identity") + 
  xlab("OE Or65a d4 vs d1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2<-ggplot(data=res_oe41[res_oe41$isexp & res_oe41$recep_or,], aes(x=symbol, y=baseMean)) +
  geom_bar(stat="identity") + 
  xlab("OE Or65a d4 vs d1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
gc <- grid.arrange(p1, p2, nrow = 1)
ggsave("newout/oe_or.pdf",gc,width = 15, height = 8)



### For OE Or65a d1 vs WT d1 -- OR
p1<-ggplot(data=res_oe_wt_d1[res_oe_wt_d1$isexp & res_oe_wt_d1$recep_or,], aes(x=symbol, y=log2FoldChange)) +
  geom_bar(stat="identity") + 
  xlab("OE Or65a d1 vs WT d1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2<-ggplot(data=res_oe_wt_d1[res_oe_wt_d1$isexp & res_oe_wt_d1$recep_or,], aes(x=symbol, y=baseMean)) +
  geom_bar(stat="identity") + 
  xlab("OE Or65a d1 vs WT d1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
gc <- grid.arrange(p1, p2, nrow = 1)
ggsave("newout/oe_wt_or_d1.pdf",gc,width = 15, height = 8)

### For OE Or65a d4 vs WT d4 -- OR
p1<-ggplot(data=res_oe_wt_d4[res_oe_wt_d4$isexp & res_oe_wt_d4$recep_or,], aes(x=symbol, y=log2FoldChange)) +
  geom_bar(stat="identity") +
  xlab("OE Or65a d4 vs WT d4") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2<-ggplot(data=res_oe_wt_d4[res_oe_wt_d4$isexp & res_oe_wt_d4$recep_or,], aes(x=symbol, y=baseMean)) +   
  geom_bar(stat="identity") +
  xlab("OE Or65a d4 vs WT d4") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
gc <- grid.arrange(p1, p2, nrow = 1)
ggsave("newout/oe_wt_or_d4.pdf",gc,width = 15, height = 8)


### For OE Or65a d4 vs d1 -- IR and GR
theylim_ir <- ylim(-1,1)
theylim_gr <- ylim(-1,1)
p1<-ggplot(data=res_oe41[res_oe41$isexp & res_oe41$recep_ir,], aes(x=symbol, y=log2FoldChange)) +
  geom_bar(stat="identity") + 
  theylim_ir +
  xlab("OE Or65a d4 vs d1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2<-ggplot(data=res_oe41[res_oe41$isexp & res_oe41$recep_gr,], aes(x=symbol, y=log2FoldChange)) +
  geom_bar(stat="identity") + 
  theylim_gr +
  xlab("OE Or65a d4 vs d1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
gc <- grid.arrange(p1, p2, nrow = 1)
ggsave("newout/oe_ir_and_gr.pdf",gc,width = 15, height = 8)



#p1<-ggplot(data=res_wt14[res_wt14$recep_or,], aes(x=symbol, y=log2FoldChange)) +
  
### For ORCO vs wt d14 -- OR
p1<-ggplot(data=res_orco_wt[res_orco_wt$isexp & res_orco_wt$recep_or & res_orco_wt$symbol!="Orco",], aes(x=symbol, y=log2FoldChange)) +
  geom_bar(stat="identity") + 
  xlab("Orco KO d14 vs WT d14") +
  ylim(-9,1.5) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p2<-ggplot(data=res_orco_wt[res_orco_wt$isexp & res_orco_wt$recep_or & res_orco_wt$symbol!="Orco",], aes(x=symbol, y=baseMean)) +
#   geom_bar(stat="identity") + 
#   xlab("Orco KO d14 vs WT d1") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2<-ggplot(data=res_orco_wt[res_orco_wt$symbol %in% intgene$symbol[intgene$category %in% c("olfactory sensory neuron markers")],], aes(x=symbol, y=log2FoldChange)) +
  geom_bar(stat="identity") + 
  xlab("Orco KO d14 vs WT d14") +
  ylim(-9,1.5) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p3<-ggplot(data=res_orco_wt[res_orco_wt$symbol %in% intgene$symbol[intgene$category %in% c("synapse markers")],], aes(x=symbol, y=log2FoldChange)) +
  geom_bar(stat="identity") + 
  xlab("Orco KO d14 vs WT d14") +
  ylim(-9,1.5) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
gc <- grid.arrange(p1, p2, p3, nrow = 1)
ggsave("newout/orco_wt_or.pdf",gc,width = 20, height = 8)




############# THE volcano

res <- merge(res_wt14,rec_count_proc, all=TRUE)
res$isexp[is.na(res$isexp)] <- TRUE
#  p1<-ggplot(data=res_wt14[res_wt14$isexp & res_wt14$recep_ir,], aes(x=symbol, y=log2FoldChange)) +
res$colcat <- "gray"
res$colcat[res$recep_ir & res$isexp] <- "red"
res$colcat[res$recep_gr & res$isexp] <- "blue"
res$colcat[res$recep_or & res$isexp] <- "green"
res$logp <- -log10(res$padj)
#res$logp <- pmin(res$logp,50)
#pdf("newout/volcano wt d1 vs d4.pdf", width = 4, height = 4)
plot(res$log2FoldChange,res$logp,col=res$colcat,pch=19, xlab="log2(Fold change)", ylab="-log10(adj.p-value)",cex=0.5, main="WT day 4 vs day 1")
res2 <- res[res$colcat!="gray",]
points(res2$log2FoldChange,res2$logp,col=res2$colcat,pch=19,cex=0.7)
#dev.off()


 ###############################


### For OE Or65a d4 vs d1 -- synapse & cilia
res <- res_oe41[res_oe41$symbol %in% intgene$symbol[intgene$category %in% c("synapse markers","synapse markers OLD","cilia")],]
p1<-ggplot(data=res, aes(x=symbol, y=log2FoldChange)) +
  geom_bar(stat="identity") + 
  xlab("OE Or65a d4 vs d1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2<-ggplot(data=res, aes(x=symbol, y=baseMean)) +
  geom_bar(stat="identity") + 
  xlab("OE Or65a d4 vs d1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
gc <- grid.arrange(p1, p2, nrow = 1)
ggsave("newout/temp.pdf",gc,width = 15, height = 8)




### For OE Or65a d4 vs d1 -- synaps & cilia
res <- res_oe41[res_oe41$symbol %in% intgene$symbol[intgene$category %in% c("tf")],]
p1<-ggplot(data=res, aes(x=symbol, y=log2FoldChange)) +
  geom_bar(stat="identity") + 
  xlab("OE Or65a d4 vs d1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2<-ggplot(data=res, aes(x=symbol, y=baseMean)) +
  geom_bar(stat="identity") + 
  xlab("OE Or65a d4 vs d1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
gc <- grid.arrange(p1, p2, nrow = 1)
ggsave("newout/temp_tf_oe.pdf",gc,width = 15, height = 8)




### For OE Or65a d4 vs d1 -- synaps & cilia
res <- res_wt14[res_wt14$symbol %in% intgene$symbol[intgene$category %in% c("tf")],]
p1<-ggplot(data=res, aes(x=symbol, y=log2FoldChange)) +
  geom_bar(stat="identity") + 
  xlab("WT d4 vs d1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2<-ggplot(data=res, aes(x=symbol, y=baseMean)) +
  geom_bar(stat="identity") + 
  xlab("WT d4 vs d1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
gc <- grid.arrange(p1, p2, nrow = 1)
ggsave("newout/temp_tf_wt.pdf",gc,width = 15, height = 8)
