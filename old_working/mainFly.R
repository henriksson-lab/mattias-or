### 19, 65,  22  --- same cell



# 
# 
# Ser lovande ut. Snabbt marscherat! Jag har spånat och de jmf som vi måste göra för lukt receptorerna är:
#   
#   
#   
#   Wt dag 1 vs 4
# Wt day 4 vs 14
# Or42b OE 1 dag vs 3 dagar
# Or42b OE 1 dag vs wt 1 dag
# 1 dag Orco mutant  vs 1 dag wt
# 1 dag Orco mutant vs 4 dagar Orco mutant


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


###################
################### Load data
###################

counts <- read.table("counts.csv", header = TRUE, stringsAsFactors = FALSE)
rownames(counts) <- counts$Geneid
counts <- counts[,-(1:6)]

cond <- read.csv("cond.csv", stringsAsFactors = FALSE, row.names = 1)
cond$treatment <- factor(cond$treatment)
#rownames(cond) == colnames(counts)
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

gene_anno$recep_ir <- str_sub(gene_anno$symbol,1,2)=="Ir" & 
  !(str_sub(gene_anno$symbol,1,3)=="Ire") & 
  !(str_sub(gene_anno$symbol,1,3)=="Irb") & 
  !(str_sub(gene_anno$symbol,1,3)=="Irk") & 
  !(str_sub(gene_anno$symbol,1,3)=="Irc") & 
  !(str_sub(gene_anno$symbol,1,3)=="Irp") & 
  !(str_sub(gene_anno$symbol,1,3)=="Iri") & 
  !(str_sub(gene_anno$symbol,1,3)=="Iru")
#sort(gene_anno$symbol[gene_anno$recep_ir])

gene_anno$recep_gr <- str_sub(gene_anno$symbol,1,2)=="Gr" & 
  !(str_sub(gene_anno$symbol,1,3)=="Gra") & 
  !(str_sub(gene_anno$symbol,1,3)=="Grd") & 
  !(str_sub(gene_anno$symbol,1,3)=="Gri") & 
  !(str_sub(gene_anno$symbol,1,3)=="Grx") 
#sort(gene_anno$symbol[gene_anno$recep_gr])

gene_anno$gclass <- ""
gene_anno$gclass[gene_anno$recep_ir] <- "Ir"
gene_anno$gclass[gene_anno$recep_or] <- "Or"
gene_anno$gclass[gene_anno$recep_gr] <- "Gr"
gene_anno$gclass[gene_anno$symbol=="Orco"] <- "Orco"
#gene_anno$gclass[gene_anno$symbol=="Or65a"] <- "Or65a"
#gene_anno$gclass[gene_anno$symbol=="Or47b"] <- "Or47b"
table(gene_anno$gclass)


geneid_oe <- gene_anno$geneid[gene_anno=="Or65a"]  #which?
#sym_receptors <- gene_anno$symbol[gene_anno$gclass=="Or"]

#geneid_receptors <- gene_anno$geneid[gene_anno$gclass=="Or"]

sym_exp_receptors_or <- read.table("expGenes.txt", stringsAsFactors = FALSE)[,1]
sym_exp_receptors_ir <- gene_anno$symbol[gene_anno$gclass=="Ir"]
sym_exp_receptors_gr <- gene_anno$symbol[gene_anno$gclass=="Gr"]

genesym_receptors <- c(sym_exp_receptors_or, sym_exp_receptors_ir, sym_exp_receptors_gr)

geneid_receptors <- gene_anno$geneid[gene_anno$symbol %in% genesym_receptors]

###################
###################


intgene <- read.csv("genelist.csv")
intgene$isin <- intgene$symbol %in% gene_anno$symbol

###################
################### Compare OE vs WT, generally
###################

keepdata <- cond$treatment!="orcdo"
dds <- DESeqDataSetFromMatrix(countData = counts[,keepdata], colData = cond[keepdata,], design = ~ treatment)

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
################### Compare octo-mutant vs WT
###################
 
keepdata <- cond$day=="14"
dds <- DESeqDataSetFromMatrix(countData = counts[,keepdata], colData = cond[keepdata,], design = ~ treatment)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast = c("treatment","orcdo","wt"))
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
################### Basal levels
###################

#rec_count <- cnt_norm[rownames(cnt_norm) %in% geneid_receptors,cond$treatment=="wt"]
rec_count <- log10(1+cnt_norm[rownames(cnt_norm) %in% geneid_receptors,cond$treatment=="wt"])
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
#plot(rec_count_proc$count14,rec_count_proc$count1)


###########3sort(unique(rec_count_proc$symbol[rec_count_proc$count4>1])) ####

p <- ggplot(rec_count_proc, aes(x = day, y = count, group = geneid, color=gclass)) + 
  geom_line() #+
  #theme(legend.position = "none") 
p

p1 <- ggplot(rec_count_proc[rec_count_proc$isinc,], aes(x = day, y = count, colour = gclass, group=geneid)) + 
  geom_line() 
  #theme(legend.position = "none") 
p2 <- ggplot(rec_count_proc[!rec_count_proc$isinc,], aes(x = day, y = count, colour = gclass, group=geneid)) + 
  geom_line() #+
#  theme(legend.position = "none") 
gc <- grid.arrange(p1, p2, nrow = 1)
ggsave("out/tc_wt.pdf",gc)

##################

p1 <- ggplot(rec_count_proc[rec_count_proc$gclass=="Ir",], aes(x = day, y = count, colour = gclass, group=geneid)) + geom_line() 
p2 <- ggplot(rec_count_proc[rec_count_proc$gclass=="Or",], aes(x = day, y = count, colour = gclass, group=geneid)) + geom_line() 
p3 <- ggplot(rec_count_proc[rec_count_proc$gclass=="Gr",], aes(x = day, y = count, colour = gclass, group=geneid)) + geom_line() 
gc <- grid.arrange(p1, p2, p3, nrow = 1)

###################


p1 <- ggplot(rec_count_proc[rec_count_proc$isinc14_4 & rec_count_proc$isinc4_1,], aes(x = day, y = count, colour = gclass, group=geneid)) + 
  geom_line() + ylim(0,5)
p2 <- ggplot(rec_count_proc[!rec_count_proc$isinc14_4 & rec_count_proc$isinc4_1,], aes(x = day, y = count, colour = gclass, group=geneid)) + 
  geom_line() + ylim(0,5)
p3 <- ggplot(rec_count_proc[rec_count_proc$isinc14_4 & !rec_count_proc$isinc4_1,], aes(x = day, y = count, colour = gclass, group=geneid)) + 
  geom_line()  + ylim(0,5)
p4 <- ggplot(rec_count_proc[!rec_count_proc$isinc14_4 & !rec_count_proc$isinc4_1,], aes(x = day, y = count, colour = gclass, group=geneid)) + 
  geom_line() + ylim(0,5)
gc <- grid.arrange(p1, p2, p3, p4, nrow = 2)
ggsave("out/tc_wt_4.pdf",gc)

recep_baseline_1  <- rec_count_proc[rec_count_proc$day==1, c("geneid","symbol","count")]
recep_baseline_14 <- rec_count_proc[rec_count_proc$day==14,c("geneid","symbol","count")]
colnames(recep_baseline_1) <- c("geneid","symbol","baseline")
colnames(recep_baseline_14) <- c("geneid","symbol","baseline")

rec_count_proc[rec_count_proc$day==14,]


# geneid        count         var day symbol organism recep
# 93  FBgn0034865 1953.7351508 346.0474942  14  Or59b     Dmel  TRUE
# 144 FBgn0038798 1632.9321069 524.4608687  14  Or92a     Dmel  TRUE
# 78  FBgn0033043 1333.6772544 262.0633342  14  Or42b     Dmel  TRUE
# 42  FBgn0026397 1296.8342942  43.2838192  14  Or22b     Dmel  TRUE
# 9   FBgn0026385  945.7815624 173.3199046  14  Or47b     Dmel  TRUE
# 63  FBgn0030204  720.4953629 182.1358130  14   Or9a     Dmel  TRUE
# 60  FBgn0030016  541.2060553  11.5712405  14   Or7a     Dmel  TRUE
# 12  FBgn0026386  506.1538940  48.0507975  14  Or47a     Dmel  TRUE
# 111 FBgn0036080  496.0730315  46.3935296  14  Or67d     Dmel  TRUE

###################
################### OE time course
###################

#rec_count <- cnt_norm[rownames(cnt_norm) %in% geneid_receptors,cond$treatment=="oe"]
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
#  theme(legend.position = "none") 
  #ylim(0, 20)
p2 <- ggplot(rec_count_proc[!rec_count_proc$isinc,], aes(x = day, y = count, colour = gclass, group=geneid)) + 
  geom_line()  + ylim(0,5)
#  theme(legend.position = "none") 
gc <- grid.arrange(p1, p2, nrow = 1)
ggsave("out/tc_oe.pdf",gc)

rec_count_proc[rec_count_proc$day==4,]

# geneid        count         var day symbol organism recep
# 114 FBgn0041625 1.444310e+04 252.3236622   4  Or65a     Dmel  TRUE  **  not in top wt   ... this is the OE gene
# 6   FBgn0026385 9.382632e+03 252.3548206   4  Or47b     Dmel  TRUE **   wt
# 96  FBgn0038798 6.465701e+02  70.1416879   4  Or92a     Dmel  TRUE
# 62  FBgn0034865 5.499769e+02  15.3322021   4  Or59b     Dmel  TRUE  wt
# 28  FBgn0026397 3.929921e+02  19.3028683   4  Or22b     Dmel  TRUE
# 52  FBgn0033043 3.318747e+02  64.9615814   4  Or42b     Dmel  TRUE  wt


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
################### barplots
###################


#wt 1 v 4
#4 vs 14
#sort by p-value? or FC. plot names


library(ggplot2)
# Basic barplot
p<-ggplot(data=df, aes(x=dose, y=len)) +
  geom_bar(stat="identity")
p

# Horizontal bar plot
p + coord_flip()






######################################################### Define expressed receptors


#sym_exp_receptors <- read.table("expGenes.txt", stringsAsFactors = FALSE)[,1]


######################################################### WT day 1 vs 4
cond$dday <- factor(sprintf("d%d",cond$day))
keepdata <- cond$treatment=="wt" & cond$day %in% c(1,4)
dds <- DESeqDataSetFromMatrix(countData = counts[,keepdata], colData = cond[keepdata,], design = ~ dday)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast = c("dday","d4","d1"))
res <- as.data.frame(res)
res$geneid <- rownames(res)
res <- merge(res, gene_anno)
res <- res[res$symbol %in% sym_exp_receptors,]
#res <- res[order(res$log2FoldChange, decreasing = FALSE),]
res <- res[order(res$baseMean, decreasing = FALSE),]

res$symbol <- factor(res$symbol, levels = res$symbol)
res_order_factor <- res$symbol

p<-ggplot(data=res, aes(x=symbol, y=log2FoldChange)) +
  geom_bar(stat="identity") + 
  coord_flip() + ylim(-2.5,2.5)
p
p_wt4_1 <- p

p<-ggplot(data=res, aes(x=symbol, y=baseMean)) +
  geom_bar(stat="identity") + 
  coord_flip() #+ ylim(-2.5,2.5)
p
p_lev <- p

######################################################### WT day 14 vs 4
cond$dday <- factor(sprintf("d%d",cond$day))
keepdata <- cond$treatment=="wt" & cond$day %in% c(4,14)
dds <- DESeqDataSetFromMatrix(countData = counts[,keepdata], colData = cond[keepdata,], design = ~ dday)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast = c("dday","d14","d4"))
res <- as.data.frame(res)
res$geneid <- rownames(res)
res <- merge(res, gene_anno)
res <- res[res$symbol %in% sym_exp_receptors,]
res <- res[order(res$log2FoldChange, decreasing = FALSE),]

res$symbol <- factor(res$symbol, levels = levels(res_order_factor))
res_order_factor <- res$symbol

p<-ggplot(data=res, aes(x=symbol, y=log2FoldChange)) +
  geom_bar(stat="identity") +
  coord_flip() + ylim(-2.5,2.5)
p_wt14_4 <- p

gc <- grid.arrange(p_wt4_1, p_wt14_4, nrow = 1)
ggsave("out/new_tc_wt.pdf",gc)




######################################################### OE day 14 vs 4
cond$dday <- factor(sprintf("d%d",cond$day))
keepdata <- cond$treatment=="oe" & cond$day %in% c(1,4)
dds <- DESeqDataSetFromMatrix(countData = counts[,keepdata], colData = cond[keepdata,], design = ~ dday)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast = c("dday","d4","d1"))
res <- as.data.frame(res)
res$geneid <- rownames(res)
res <- merge(res, gene_anno)
res <- res[res$symbol %in% sym_exp_receptors,]
res <- res[order(res$log2FoldChange, decreasing = FALSE),]

res$symbol <- factor(res$symbol, levels = levels(res_order_factor))
res_order_factor <- res$symbol

p<-ggplot(data=res, aes(x=symbol, y=log2FoldChange)) +
  geom_bar(stat="identity") +
  coord_flip() + ylim(-2.5,2.5)
p_oe <- p
p_oe
#gc <- grid.arrange(p_wt4_1, p_wt14_4, p_oe, nrow = 1)
gc <- grid.arrange(p_lev, p_wt4_1, p_wt14_4, p_oe, nrow = 1)
ggsave("out/new_tc_oe.pdf",gc)

