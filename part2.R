### low-level function
plotxy.2 <- function(colX,colY, titColX, titColY, xlim=c(0,5), ylim=c(0,5), title=""){
  plot(rec_count[,colX], rec_count[,colY], pch=19, col=rec_count$colcat, cex=0.5, 
       xlim=xlim,ylim=ylim,
       xlab=titColX,
       ylab=titColY,
       main=title)
  temp <- rec_count[rec_count$colcat!="gray",]
  lines(c(-7,7),c(-7,7))
  points(temp[,colX], temp[,colY], pch=19, col=temp$colcat)
}

### make scatter plots
plotxy <- function(colX,colY, titColX, titColY, fname,xlim=c(0,5), ylim=c(0,5), title=""){
  plotxy.2(colX,colY,titColX, titColY, xlim=xlim, ylim=ylim, title=title)
  pdf(fname, w=5, h=5)
  plotxy.2(colX,colY,titColX, titColY, xlim=xlim, ylim=ylim, title=title)
  dev.off()
  pdf(sprintf("%s_noborder.pdf",fname), w=5, h=5)
  plotxy.2(colX,colY,"", "", xlim=xlim, ylim=ylim, title="")
  dev.off()
}

### scatter flow, reverse x and y
plotxy.rev <- function(colX,colY, titColX, titColY, fname,xlim=c(0,5), ylim=c(0,5), title=""){
  plotxy(colY,colX, titColY, titColX, fname,xlim=xlim, ylim=ylim, title=title)
}





###################
################### Plots of wildtype vs orco
###################

if(TRUE){
  #SF normalized
  rec_count <- cnt_norm
  rec_count <- log10(1+rec_count)
  
} else {
  #CPM normalized
  rec_count <- cnt_norm
  for(i in 1:ncol(rec_count)){
    rec_count[,i] <- rec_count[,i]/(sum(rec_count[,i])/1e6)
  }
  rec_count <- log10(1+rec_count)
}


#Merge conditions
newcond <- unique(cond[,c("day","treatment")])
newcounts <- matrix(nrow=nrow(rec_count), ncol=nrow(newcond))
for(i in 1:ncol(newcounts)){
  newcounts[,i] <- rowMeans(rec_count[,cond$day==newcond$day[i] & cond$treatment==newcond$treatment[i], drop=FALSE])
}
colnames(newcounts) <- sprintf("%s_%s",newcond$treatment,newcond$day)
rownames(newcond) <- colnames(newcounts)
rownames(newcounts) <- rownames(counts)


rec_count <- as.data.frame(newcounts)
rec_count$geneid <- rownames(rec_count)
rec_count$colcat <- "gray"
rec_count$colcat[rec_count$geneid %in% gene_anno[gene_anno$is_trichoids,]$geneid] <- "darkmagenta"
rec_count$colcat[rec_count$geneid %in% gene_anno[gene_anno$is_basiconic,]$geneid] <- "darkgreen"


rec_count$delta_orco14 <- rec_count$orco_4 - rec_count$orco_1
rec_count$delta_wt14 <- rec_count$wt_4 - rec_count$wt_1



# 
# ################
# plotxy.rev("wt_1","wt_4","Log10(WT day 1 CPM)","Log10(WT day 4 CPM)", "orco_out/1_1_wt4_vs_wt1.pdf", title = "Control 1-4 days ORs")
# plotxy.rev("oe_1","oe_4","Log10(OE day 1 CPM)","Log10(OE day 4 CPM)", "orco_out/1_2_oe4_vs_oe1.pdf", title = "OE 1-4 days ORs")
# plotxy.rev("oe_1","wt_1","Log10(OE day 1 CPM)","Log10(WT day 1 CPM)", "orco_out/1_3_oe1_vs_wt1.pdf", title = "OE vs Control 1 day ORs")
# 
# plotxy.rev("oe_4","wt_4","Log10(OE day 4 CPM)","Log10(wt day 4 CPM)", "orco_out/2_1_oe4_vs_wt4.pdf", title = "OE vs Control 4 day ORs")
# plotxy.rev("orco_1","orco_4","Log10(ORCO day 1 CPM)","Log10(ORCO day 4 CPM)", "orco_out/2_2_orco1_vs_orco4.pdf", title = "Orco 1-4 days ORs")
# plotxy.rev("orco_1","wt_1","Log10(ORCO day 1 CPM)","Log10(WT day 1 CPM)", "orco_out/2_3_orco1_vs_orco4.pdf", title = "Orco vs control 1 day ORs")
# 
# plotxy.rev("orco_4","wt_4","Log10(ORCO day 4 CPM)","Log10(WT day 4 CPM)", "orco_out/3_1_orco4_vs_wt4.pdf", title = "Orco vs control 4 day ORs")

###########################################################
############# Expressed genes? ############################
###########################################################


#### Decide which genes are "expressed"
cutoff_gr <- rec_count$wt_4[rec_count$geneid %in% gene_anno[gene_anno$symbol=="Gr28a",]$geneid]
cutoff_ir <- rec_count$wt_4[rec_count$geneid %in% gene_anno[gene_anno$symbol=="Ir60a",]$geneid]
cutoff_or <- rec_count$wt_4[rec_count$geneid %in% gene_anno[gene_anno$symbol=="Or65c",]$geneid]
cutoff_general <- 1

rec_count$is_exp <- TRUE
rec_count$is_exp[rec_count$geneid %in% gene_anno[gene_anno$recep_gr,]$geneid & rec_count$wt_4<cutoff_gr] <- FALSE
rec_count$is_exp[rec_count$geneid %in% gene_anno[gene_anno$recep_ir,]$geneid & rec_count$wt_4<cutoff_ir] <- FALSE
rec_count$is_exp[rec_count$geneid %in% gene_anno[gene_anno$recep_or,]$geneid & rec_count$wt_4<cutoff_or] <- FALSE

rec_count$is_exp[rec_count$geneid %in% gene_anno[gene_anno$other_int_obp,]$geneid & rec_count$wt_4<cutoff_general] <- FALSE

rec_count$is_exp[rec_count$geneid %in% gene_anno[gene_anno$symbol=="Or19a",]$geneid] <- TRUE




###
#write.csv(merge(gene_anno,rec_count),"out/de_table.csv")

###########################################################
############# DE tests ####################################
###########################################################


##############################
keep <- cond$day==1 & cond$treatment %in% c("wt","orco")
dds <- DESeqDataSetFromMatrix(countData = counts[,keep], colData = cond[keep,], design = ~treatment)
dds <- DESeq(dds)
r <- as.data.frame(results(dds, contrast = c("treatment","orco","wt")))
r$geneid <- rownames(r)
r <- merge(r, unique(gene_anno[,c("symbol","geneid")]))
r <- r[order(r$pvalue),]
write.csv(r,"out/de/de_day1_wt_orco.csv")
##############################
keep <- cond$day==4 & cond$treatment %in% c("wt","orco")
dds <- DESeqDataSetFromMatrix(countData = counts[,keep], colData = cond[keep,], design = ~treatment)
dds <- DESeq(dds)
r <- as.data.frame(results(dds, contrast = c("treatment","orco","wt")))
r$geneid <- rownames(r)
r <- merge(r, unique(gene_anno[,c("symbol","geneid")]))
r <- r[order(r$pvalue),]
write.csv(r,"out/de/de_day4_wt_orco.csv")
##############################
keep <- cond$day==14 & cond$treatment %in% c("wt","orco")
dds <- DESeqDataSetFromMatrix(countData = counts[,keep], colData = cond[keep,], design = ~treatment)
dds <- DESeq(dds)
r <- as.data.frame(results(dds, contrast = c("treatment","orco","wt")))
r$geneid <- rownames(r)
r <- merge(r, unique(gene_anno[,c("symbol","geneid")]))
r <- r[order(r$pvalue),]
write.csv(r,"out/de/de_day14_wt_orco.csv")






##############################
keep <- cond$day==1 & cond$treatment %in% c("wt","oe")
dds <- DESeqDataSetFromMatrix(countData = counts[,keep], colData = cond[keep,], design = ~treatment)
dds <- DESeq(dds)
r <- as.data.frame(results(dds, contrast = c("treatment","oe","wt")))
r$geneid <- rownames(r)
r <- merge(r, unique(gene_anno[,c("symbol","geneid")]))
r <- r[order(r$pvalue),]
write.csv(r,"out/de/de_day1_wt_oe.csv")
##############################
keep <- cond$day==4 & cond$treatment %in% c("wt","oe")
dds <- DESeqDataSetFromMatrix(countData = counts[,keep], colData = cond[keep,], design = ~treatment)
dds <- DESeq(dds)
r <- as.data.frame(results(dds, contrast = c("treatment","oe","wt")))
r$geneid <- rownames(r)
r <- merge(r, unique(gene_anno[,c("symbol","geneid")]))
r <- r[order(r$pvalue),]
write.csv(r,"out/de/de_day4_wt_oe.csv")
##############################
keep <- cond$day==14 & cond$treatment %in% c("wt","oe")
dds <- DESeqDataSetFromMatrix(countData = counts[,keep], colData = cond[keep,], design = ~treatment)
dds <- DESeq(dds)
r <- as.data.frame(results(dds, contrast = c("treatment","oe","wt")))
r$geneid <- rownames(r)
r <- merge(r, unique(gene_anno[,c("symbol","geneid")]))
r <- r[order(r$pvalue),]
write.csv(r,"out/de/de_day14_wt_oe.csv")



##############################
keep <- cond$day %in% c(1,4) & cond$treatment %in% c("wt")
cond$dday <- sprintf("d%s",cond$day)
dds <- DESeqDataSetFromMatrix(countData = counts[,keep], colData = cond[keep,], design = ~dday) 
dds <- DESeq(dds)
r <- as.data.frame(results(dds))
r <- as.data.frame(results(dds, contrast = c("dday","d4","d1")))
r$geneid <- rownames(r)
r <- merge(r, unique(gene_anno[,c("symbol","geneid")]))
r <- r[order(r$pvalue),]
write.csv(r,"out/de/de_wt day4 vs day1.csv")  
##############################
keep <- cond$day %in% c(4,14) & cond$treatment %in% c("wt")
cond$dday <- sprintf("d%s",cond$day)
dds <- DESeqDataSetFromMatrix(countData = counts[,keep], colData = cond[keep,], design = ~dday) 
dds <- DESeq(dds)
r <- as.data.frame(results(dds))
r <- as.data.frame(results(dds, contrast = c("dday","d14","d4")))
r$geneid <- rownames(r)
r <- merge(r, unique(gene_anno[,c("symbol","geneid")]))
r <- r[order(r$pvalue),]
write.csv(r,"out/de/de_wt day14 vs day4.csv")  


##############################
keep <- cond$day %in% c(1,4) & cond$treatment %in% c("orco")
cond$dday <- sprintf("d%s",cond$day)
dds <- DESeqDataSetFromMatrix(countData = counts[,keep], colData = cond[keep,], design = ~dday) 
dds <- DESeq(dds)
r <- as.data.frame(results(dds))
r <- as.data.frame(results(dds, contrast = c("dday","d4","d1")))
r$geneid <- rownames(r)
r <- merge(r, unique(gene_anno[,c("symbol","geneid")]))
r <- r[order(r$pvalue),]
write.csv(r,"out/de/de_orco day4 vs day1.csv")  
##############################
keep <- cond$day %in% c(4,14) & cond$treatment %in% c("orco")
cond$dday <- sprintf("d%s",cond$day)
dds <- DESeqDataSetFromMatrix(countData = counts[,keep], colData = cond[keep,], design = ~dday) 
dds <- DESeq(dds)
r <- as.data.frame(results(dds))
r <- as.data.frame(results(dds, contrast = c("dday","d14","d4")))
r$geneid <- rownames(r)
r <- merge(r, unique(gene_anno[,c("symbol","geneid")]))
r <- r[order(r$pvalue),]
write.csv(r,"out/de/de_orco day14 vs day4.csv")  





##############################
keep <- cond$day %in% c(1,4) & cond$treatment %in% c("oe")
cond$dday <- sprintf("d%s",cond$day)
dds <- DESeqDataSetFromMatrix(countData = counts[,keep], colData = cond[keep,], design = ~dday) 
dds <- DESeq(dds)
r <- as.data.frame(results(dds))
r <- as.data.frame(results(dds, contrast = c("dday","d4","d1")))
r$geneid <- rownames(r)
r <- merge(r, unique(gene_anno[,c("symbol","geneid")]))
r <- r[order(r$pvalue),]
write.csv(r,"out/de/de_oe day4 vs day1.csv")  

##############################
keep <- cond$day %in% c(4,14) & cond$treatment %in% c("oe")
cond$dday <- sprintf("d%s",cond$day)
dds <- DESeqDataSetFromMatrix(countData = counts[,keep], colData = cond[keep,], design = ~dday) 
dds <- DESeq(dds)
r <- as.data.frame(results(dds))
r <- as.data.frame(results(dds, contrast = c("dday","d14","d4")))
r$geneid <- rownames(r)
r <- merge(r, unique(gene_anno[,c("symbol","geneid")]))
r <- r[order(r$pvalue),]
write.csv(r,"out/de/de_oe day14 vs day4.csv")  





# ##############################
# keep <- cond$day %in% c(1,4) & cond$treatment %in% c("orco","wt")
# cond$dday_treatment <- sprintf("d%s_%s",cond$day, cond$treatment)
# dds <- DESeqDataSetFromMatrix(countData = counts[,keep], colData = cond[keep,], design = ~dday_treatment+0)
# dds <- DESeq(dds)
# #resultsNames(dds) #"dday_treatmentd1_orco" "dday_treatmentd1_wt"   "dday_treatmentd4_orco" "dday_treatmentd4_wt"
# r <- as.data.frame(results(dds, contrast = c(1,-1,-1,1)))
# r$geneid <- rownames(r)
# r <- merge(r, unique(gene_anno[,c("symbol","geneid")]))
# r <- r[order(r$pvalue),]
# write.csv(r,"out/de/de_delta.csv")









###########################################################
############# Figure 1 ####################################
###########################################################

plot_rec_volcano <- function(fname, title, outfname){
  pdf(outfname, w=5, h=5)
  
  x<-read.csv("out/de/de_wt day4 vs day1.csv", stringsAsFactors = FALSE)  
  x<-merge(x,data.frame(
    geneid=rownames(rec_count),
    colcat=rec_count$colcat))
  plot(x$log2FoldChange, -log10(x$padj), pch=19, col=as.character(x$colcat), cex=0.5, 
       xlab="Log2 fold change",
       ylab="-Log10 p.adj",
       main=title)
  x <- x[x$colcat!="gray",]
  points(x$log2FoldChange, -log10(x$padj), pch=19, col=as.character(x$colcat))
  dev.off()
  
  sum(x$colcat!="gray" & x$padj<0.01)
}


## only the wt cells; 1-4  4-14  DPE
#b: ORs GRs IRs
#c: synapse markers
#d: cilia transport (IFT and BBS)
#e: OBPs

############## b
rec_count$colcat <- "gray"
rec_count$colcat[rec_count$is_exp & rec_count$geneid %in% gene_anno[gene_anno$recep_ir,]$geneid] <- "red"
rec_count$colcat[rec_count$is_exp & rec_count$geneid %in% gene_anno[gene_anno$recep_gr,]$geneid] <- "blue"
rec_count$colcat[rec_count$is_exp & rec_count$geneid %in% gene_anno[gene_anno$recep_or,]$geneid] <- "darkgreen"

plotxy("wt_1","wt_4", "wt 1 DPE","wt 4 DPE",  "out/fig1/b_wt4_vs_wt1_ORsGRsIRs.pdf",  title = "1 DPE vs 4 DPE")
plotxy("wt_4","wt_14","wt 4 DPE","wt 14 DPE", "out/fig1/b_wt14_vs_wt4_ORsGRsIRs.pdf", title = "4 DPE vs 14 DPE")


plot_rec_volcano("out/de/de_wt day4 vs day1.csv", "1 DPE vs 4 DPE", "out/volcano/fig1_b_wt4_vs_wt1_ORsGRsIRs.pdf")   #51
plot_rec_volcano("out/de/de_wt day14 vs day4.csv", "4 DPE vs 14 DPE", "out/volcano/fig1_b_wt14_vs_wt4_ORsGRsIRs.pdf") #51


############## c
rec_count$colcat <- "gray"
rec_count$colcat[rec_count$is_exp & rec_count$geneid %in% gene_anno[gene_anno$is_synapse,]$geneid] <- "red"
#rec_count$colcat[rec_count$is_exp & rec_count$geneid %in% gene_anno[gene_anno$symbol=="elav",]$geneid] <- "blue"


plotxy("wt_1","wt_4", "wt 1 DPE","wt 4 DPE",  "out/fig1/c_wt4_vs_wt1_synapse.pdf",  title = "1 DPE vs 4 DPE")
plotxy("wt_4","wt_14","wt 4 DPE","wt 14 DPE", "out/fig1/c_wt14_vs_wt4_synapse.pdf", title = "4 DPE vs 14 DPE")

plot_rec_volcano("out/de/de_wt day4 vs day1.csv", "1 DPE vs 4 DPE", "out/volcano/fig1_c_wt4_vs_wt1_synapse.pdf") #10
plot_rec_volcano("out/de/de_wt day14 vs day4.csv", "4 DPE vs 14 DPE", "out/volcano/fig1_c_wt14_vs_wt4_synapse.pdf") #10


############## d
rec_count$colcat <- "gray"
rec_count$colcat[rec_count$is_exp & rec_count$geneid %in% gene_anno[gene_anno$is_cilia,]$geneid] <- "red"

plotxy("wt_1","wt_4", "wt 1 DPE","wt 4 DPE",  "out/fig1/d_wt4_vs_wt1_cilia.pdf",  title = "1 DPE vs 4 DPE")
plotxy("wt_4","wt_14","wt 4 DPE","wt 14 DPE", "out/fig1/d_wt14_vs_wt4_cilia.pdf", title = "4 DPE vs 14 DPE")

plot_rec_volcano("out/de/de_wt day4 vs day1.csv", "1 DPE vs 4 DPE", "out/volcano/fig1_d_wt4_vs_wt1_cilia.pdf")  #13
plot_rec_volcano("out/de/de_wt day14 vs day4.csv", "4 DPE vs 14 DPE", "out/volcano/fig1_d_wt14_vs_wt4_cilia.pdf") #13

############## e
rec_count$colcat <- "gray"
rec_count$colcat[rec_count$is_exp & rec_count$geneid %in% gene_anno[gene_anno$other_int_obp,]$geneid] <- "red"

plotxy("wt_1","wt_4", "wt 1 DPE","wt 4 DPE",  "out/fig1/e_wt4_vs_wt1_OBP.pdf",  title = "1 DPE vs 4 DPE")
plotxy("wt_4","wt_14","wt 4 DPE","wt 14 DPE", "out/fig1/e_wt14_vs_wt4_OBP.pdf", title = "4 DPE vs 14 DPE")

plot_rec_volcano("out/de/de_wt day4 vs day1.csv", "1 DPE vs 4 DPE", "out/volcano/fig1_e_wt4_vs_wt1_OBP.pdf") #23
plot_rec_volcano("out/de/de_wt day14 vs day4.csv", "4 DPE vs 14 DPE", "out/volcano/fig1_e_wt4_vs_wt1_OBP.pdf") #23



###########################################################
############# Figure 3 ####################################
###########################################################

rec_count$colcat <- "gray"
rec_count$colcat[rec_count$geneid %in% gene_anno[gene_anno$is_trichoids,]$geneid] <- "darkmagenta"
rec_count$colcat[rec_count$geneid %in% gene_anno[gene_anno$is_basiconic,]$geneid] <- "darkgreen"

## control 1-4    OE 1-4    Orco- 1-4
#y: day4
#x: day1

plotxy("wt_1","wt_4","Log10(WT day 1 CPM)","Log10(WT day 4 CPM)", "out/fig3/wt4_vs_wt1.pdf", title = "Control 1-4 days ORs")
plotxy("oe_1","oe_4","Log10(OE day 1 CPM)","Log10(OE day 4 CPM)", "out/fig3/oe4_vs_oe1.pdf", title = "OE 1-4 days ORs")
#plotxy.rev("oe_1","wt_1","Log10(OE day 1 CPM)","Log10(WT day 1 CPM)", "orco_out/1_3_oe1_vs_wt1.pdf", title = "OE vs Control 1 day ORs")
plotxy("orco_1","orco_4","Log10(ORCO day 1 CPM)","Log10(ORCO day 4 CPM)", "out/fig3/orco1_vs_orco4.pdf", title = "Orco 1-4 days ORs")


plot_rec_volcano("out/de/de_wt day4 vs day1.csv", "Control 1-4 days ORs", "out/volcano/fig3_wt4_vs_wt1.pdf")  #21
plot_rec_volcano("out/de/de_oe day4 vs day1.csv", "OE 1-4 days ORs", "out/volcano/fig3_oe4_vs_oe1.pdf") #21
plot_rec_volcano("out/de/de_orco day4 vs day1.csv", "Orco 1-4 days ORs", "out/volcano/fig3_orco1_vs_orco4.pdf") #21


###########################################################
############# Sup Figure 1 ################################
###########################################################

rec_count$colcat <- "gray"
rec_count$colcat[rec_count$geneid %in% gene_anno[gene_anno$is_trichoids,]$geneid] <- "darkmagenta"
rec_count$colcat[rec_count$geneid %in% gene_anno[gene_anno$is_basiconic,]$geneid] <- "darkgreen"

####### wt d1    wt d4    ????????????
#OE
#ORco-

plotxy.rev("oe_1","wt_4","Peb-Gal4; UAS-Or42b","WT 1 DPE", "out/sup1/oe1_vs_wt1.pdf", title = "OE vs Control 4 day ORs")
plotxy.rev("oe_4","wt_4","Peb-Gal4; UAS-Or42b","WT 4 DPE", "out/sup1/oe4_vs_wt4.pdf", title = "OE vs Control 4 day ORs")
plotxy.rev("orco_1","wt_1","Orco-/-","WT 1 DPE", "out/sup1/orco1_vs_wt1.pdf", title = "Orco vs control 1 day ORs")
plotxy.rev("orco_4","wt_4","Orco-/-","WT 4 DPE", "out/sup1/orco4_vs_wt4.pdf", title = "Orco vs control 4 day ORs")  



###########################################################
###########################################################
###########################################################
###########################################################



rec_count$colcat <- "gray"
rec_count$colcat[rec_count$is_exp & rec_count$geneid %in% gene_anno[gene_anno$is_synapse,]$geneid] <- "red"

plotxy.rev("orco_4","wt_4","Orco-/-","WT 4 DPE", "out/sup1/orco4_vs_wt4_synapse.pdf", title = "Orco vs control 4 day synapse")



###########################################################
############# Summary CSV file ############################
###########################################################

################### First summary file

tosave <- rec_count
tosave <- merge(tosave,gene_anno[,c("symbol", "geneid")])
tosave <- tosave[,c("geneid","symbol",
  "wt_1","orco_1","oe_1",
  "wt_4","orco_4","oe_4")]
colnames(tosave) <- c("geneid","symbol",
  "Level Control DPE1","Level Orco DPE1","Level Or42b DPE1",
  "Level Control DPE4","Level Orco DPE4","Level Or42b DPE4"
  )
tosave[,str_starts(colnames(tosave),"Level")] <- signif(tosave[,str_starts(colnames(tosave),"Level")], digits=3)

read_and_rename <- function(fname,thecol){
  d <- read.csv(fname, stringsAsFactors = FALSE)
  dd <- data.frame(geneid=d$geneid)#matrix(ncol=3, nrow=nrow(d))
  dd[,paste("FC",thecol)] <-  signif(d$log2FoldChange,digits=3)
  dd[,paste("pvalue",thecol)]  <- signif(d$pvalue,digits=3)
  dd
}

tosave_comparison <- merge(
    read_and_rename("out/de/de_day1_wt_orco.csv","Orco/Control DPE1"),
    read_and_rename("out/de/de_day4_wt_orco.csv","Orco/Control DPE4"))
tosave_comparison <- merge(
  read_and_rename("out/de/de_day1_wt_oe.csv","Or42b/Control DPE1"),
  read_and_rename("out/de/de_day4_wt_oe.csv","Or42b/Control DPE4"))

tosave_comparison <- merge(tosave_comparison,read_and_rename("out/de/de_orco day4 vs day1.csv","Orco DPE4/DPE1"))
#tosave_comparison <- merge(tosave_comparison,read_and_rename("out/de/de_orco day14 vs day4.csv","Orco DPE14/DPE4"))
tosave_comparison <- merge(tosave_comparison,read_and_rename("out/de/de_oe day4 vs day1.csv","Or42b DPE4/DPE1"))
#tosave_comparison <- merge(tosave_comparison,read_and_rename("out/de/de_oe day14 vs day4.csv","Or42b DPE14/DPE4"))


tosave_comparison <- merge(tosave, tosave_comparison)

## only interesting genes
genes_to_save <- gene_anno[gene_anno$recep_or,]$geneid
tosave_comparison_red <- tosave_comparison[tosave_comparison$geneid %in% genes_to_save,]

## exclude non-exp genes  
exp_geneid <- rec_count$geneid[rec_count$is_exp]
tosave_comparison_red <- tosave_comparison_red[tosave_comparison_red$geneid %in% exp_geneid,]
#genes_to_save <- intersect(genes_to_save,exp_geneid)


#tosave_comparison_or <- tosave_comparison[tosave_comparison$geneid %in% gene_anno[gene_anno$recep_or,]$geneid,]
#tosave_comparison_or <- na.omit(tosave_comparison_or)
tosave_comparison_red <- tosave_comparison_red[order(tosave_comparison_red$symbol),]

write.csv(tosave_comparison_red, "out/summary data PERTURB red.csv",row.names = FALSE)



############################# second summary file 

tosave <- rec_count
tosave <- merge(tosave,gene_anno[,c("symbol", "geneid")])
tosave <- tosave[,c("geneid","symbol",
                    "wt_1","wt_4","wt_14")]
colnames(tosave) <- c("geneid","symbol",
                      "Level Control DPE1","Level Control DPE4","Level Control DPE14"
)
tosave[,str_starts(colnames(tosave),"Level")] <- signif(tosave[,str_starts(colnames(tosave),"Level")], digits=3)

tosave_comparison <- merge(read_and_rename("out/de/de_wt day4 vs day1.csv","Control DPE4/DPE1"),
                           read_and_rename("out/de/de_wt day14 vs day4.csv","Control DPE14/DPE4"))
tosave_comparison <- merge(tosave, tosave_comparison)

tosave_comparison$genecat <- ""
tosave_comparison$genecat[tosave_comparison$geneid %in% gene_anno[gene_anno$recep_or,]$geneid] <- "Receptor"
tosave_comparison$genecat[tosave_comparison$geneid %in% gene_anno[gene_anno$recep_gr,]$geneid] <- "Receptor"
tosave_comparison$genecat[tosave_comparison$geneid %in% gene_anno[gene_anno$recep_ir,]$geneid] <- "Receptor"
tosave_comparison$genecat[tosave_comparison$geneid %in% gene_anno[gene_anno$other_int_obp,]$geneid] <- "OBP"
tosave_comparison$genecat[tosave_comparison$geneid %in% gene_anno[gene_anno$is_synapse,]$geneid] <- "Synapse"

# genes_to_save <- c(
#   gene_anno[gene_anno$recep_or,]$geneid, 
#   gene_anno[gene_anno$recep_gr,]$geneid,
#   gene_anno[gene_anno$recep_ir,]$geneid,
#   gene_anno[gene_anno$other_int_obp,]$geneid,
#   gene_anno[gene_anno$is_synapse_osn,]$geneid)
# tosave_comparison_red <- tosave_comparison[tosave_comparison$geneid %in% genes_to_save,]
tosave_comparison_red <- tosave_comparison[tosave_comparison$genecat!="",]

## only interesting genes
genes_to_save <- gene_anno[gene_anno$recep_or,]$geneid

## exclude non-exp genes  
exp_geneid <- rec_count$geneid[rec_count$is_exp]
tosave_comparison_red <- tosave_comparison_red[tosave_comparison_red$geneid %in% exp_geneid,]

#tosave_comparison_red <- tosave_comparison[tosave_comparison$geneid %in% genes_to_save,]
#tosave_comparison_or <- na.omit(tosave_comparison_or)
tosave_comparison_red <- tosave_comparison_red[order(tosave_comparison_red$symbol),]
tosave_comparison_red <- tosave_comparison_red[order(tosave_comparison_red$genecat),]

write.csv(tosave_comparison_red, "out/summary data CONTROL red.csv",row.names = FALSE)
