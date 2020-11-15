plotxy.2 <- function(colX,colY, titColX, titColY, xlim=c(0,5), ylim=c(0,5), title=""){
  # rec_count_proc$colcat <- "gray"
  # rec_count_proc$colcat[rec_count_proc$recep_ir & rec_count_proc$isexp] <- "red"
  # rec_count_proc$colcat[rec_count_proc$recep_gr & rec_count_proc$isexp] <- "blue"
  # rec_count_proc$colcat[rec_count_proc$recep_or & rec_count_proc$isexp] <- "green"
  # rec_count_proc$colcat[rec_count_proc$recep_or_orig] <- "green"
  
  plot(rec_count_proc[,colX], rec_count_proc[,colY], pch=19, col=rec_count_proc$colcat, cex=0.5, 
       xlim=xlim,ylim=ylim,
       xlab=titColX,
       ylab=titColY,
       main=title)
  temp <- rec_count_proc[rec_count_proc$colcat!="gray",]
  lines(c(-7,7),c(-7,7))
  points(temp[,colX], temp[,colY], pch=19, col=temp$colcat)
}

plotxy <- function(colX,colY, titColX, titColY, fname,xlim=c(0,5), ylim=c(0,5), title=""){
  plotxy.2(colX,colY,titColX, titColY, xlim=xlim, ylim=ylim, title=title)
  pdf(fname, w=5, h=5)
  plotxy.2(colX,colY,titColX, titColY, xlim=xlim, ylim=ylim, title=title)
  dev.off()
}

plotxy.rev <- function(colX,colY, titColX, titColY, fname,xlim=c(0,5), ylim=c(0,5), title=""){
  plotxy(colY,colX, titColY, titColX, fname,xlim=xlim, ylim=ylim, title=title)
}





###################
################### Plots of wildtype vs orco
###################

#deseq normalized
#rec_count <- log10(1+cnt_norm)

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
#rec_count$colcat[rec_count$geneid %in% gene_anno[gene_anno$recep_or,]$geneid] <- "darkgreen"

rec_count$colcat[rec_count$geneid %in% gene_anno[gene_anno$is_trichoids,]$geneid] <- "magenta"
rec_count$colcat[rec_count$geneid %in% gene_anno[gene_anno$is_basiconic,]$geneid] <- "darkgreen"

# gene_anno$is_trichoids <- gene_anno$symbol %in% intgene$symbol[intgene$category %in% c("trichoids")]
# gene_anno$is_basiconic <- gene_anno$symbol %in% intgene$symbol[intgene$category %in% c("basiconic")]



rec_count$delta_orco14 <- rec_count$orco_4-rec_count$orco_1
rec_count$delta_wt14 <- rec_count$wt_4-rec_count$wt_1

rec_count_proc <- rec_count

################
plotxy.rev("wt_1","wt_4","Log10(WT day 1 CPM)","Log10(WT day 4 CPM)", "orco_out/1_1_wt4_vs_wt1.pdf", title = "Control 1-4 days ORs")
plotxy.rev("oe_1","oe_4","Log10(OE day 1 CPM)","Log10(OE day 4 CPM)", "orco_out/1_2_oe4_vs_oe1.pdf", title = "OE 1-4 days ORs")
plotxy.rev("oe_1","wt_1","Log10(OE day 1 CPM)","Log10(WT day 1 CPM)", "orco_out/1_3_oe1_vs_wt1.pdf", title = "OE vs Control 1 day ORs")

plotxy.rev("oe_4","wt_4","Log10(OE day 4 CPM)","Log10(wt day 4 CPM)", "orco_out/2_1_oe4_vs_wt4.pdf", title = "OE vs Control 4 day ORs")
plotxy.rev("orco_1","orco_4","Log10(ORCO day 1 CPM)","Log10(ORCO day 4 CPM)", "orco_out/2_2_orco1_vs_orco4.pdf", title = "Orco 1-4 days ORs")
plotxy.rev("orco_1","wt_1","Log10(ORCO day 1 CPM)","Log10(WT day 1 CPM)", "orco_out/2_3_orco1_vs_orco4.pdf", title = "Orco vs control 1 day ORs")

plotxy.rev("orco_4","wt_4","Log10(ORCO day 4 CPM)","Log10(WT day 4 CPM)", "orco_out/3_1_orco4_vs_wt4.pdf", title = "Orco vs control 4 day ORs")






###
write.csv(merge(gene_anno,rec_count),"orco_out/de_table.csv")


##############################
keep <- cond$day==1 & cond$treatment %in% c("wt","orco")
dds <- DESeqDataSetFromMatrix(countData = counts[,keep], colData = cond[keep,], design = ~treatment)
dds <- DESeq(dds)
r <- as.data.frame(results(dds))
r$geneid <- rownames(r)
r <- merge(r, unique(gene_anno[,c("symbol","geneid")]))
r <- r[order(r$pvalue),]
#t(log10(1+cnt_norm)["FBgn0037324",keep,drop=FALSE])
write.csv(r,"orco_out/de_day1_wt_orco.csv")
##############################
keep <- cond$day==4 & cond$treatment %in% c("wt","orco")
dds <- DESeqDataSetFromMatrix(countData = counts[,keep], colData = cond[keep,], design = ~treatment)
dds <- DESeq(dds)
r <- as.data.frame(results(dds, contrast = c("treatment","wt","orco")))
r$geneid <- rownames(r)
r <- merge(r, unique(gene_anno[,c("symbol","geneid")]))
r <- r[order(r$pvalue),]
#t(log10(1+cnt_norm)["FBgn0037324",keep,drop=FALSE])
write.csv(r,"orco_out/de_day4_wt_orco.csv")
##############################
keep <- cond$day==14 & cond$treatment %in% c("wt","orco")
dds <- DESeqDataSetFromMatrix(countData = counts[,keep], colData = cond[keep,], design = ~treatment)
dds <- DESeq(dds)
r <- as.data.frame(results(dds, contrast = c("treatment","wt","orco")))
r$geneid <- rownames(r)
r <- merge(r, unique(gene_anno[,c("symbol","geneid")]))
r <- r[order(r$pvalue),]
t(log10(1+cnt_norm)["FBgn0037324",keep,drop=FALSE])
write.csv(r,"orco_out/de_day14_wt_orco.csv")
##############################
keep <- cond$day %in% c(1,4) & cond$treatment %in% c("orco")
cond$dday <- sprintf("d%s",cond$day)
dds <- DESeqDataSetFromMatrix(countData = counts[,keep], colData = cond[keep,], design = ~dday)  ############### beeep. wrong!
dds <- DESeq(dds)
r <- as.data.frame(results(dds))
r$geneid <- rownames(r)
r <- merge(r, unique(gene_anno[,c("symbol","geneid")]))
r <- r[order(r$pvalue),]
write.csv(r,"orco_out/de_day4orco_day1orco.csv")   #was wrong!












##############################
keep <- cond$day %in% c(1,4) & cond$treatment %in% c("orco","wt")
cond$dday_treatment <- sprintf("d%s_%s",cond$day, cond$treatment)
dds <- DESeqDataSetFromMatrix(countData = counts[,keep], colData = cond[keep,], design = ~dday_treatment+0)
dds <- DESeq(dds)
#resultsNames(dds) #"dday_treatmentd1_orco" "dday_treatmentd1_wt"   "dday_treatmentd4_orco" "dday_treatmentd4_wt"
r <- as.data.frame(results(dds, contrast = c(1,-1,-1,1)))
r$geneid <- rownames(r)
r <- merge(r, unique(gene_anno[,c("symbol","geneid")]))
r <- r[order(r$pvalue),]
write.csv(r,"orco_out/de_delta.csv")


