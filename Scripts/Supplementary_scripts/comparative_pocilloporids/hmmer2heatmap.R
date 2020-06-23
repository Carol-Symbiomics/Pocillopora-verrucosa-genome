
library(rhmmer)
library(gplots)
library(matrixStats)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/Pver_genome/")

##############################
# Process HMMsearch results #
##############################

pdam=read_domtblout("Pdam_domains")
pdam.f=subset(pdam, sequence_evalue < 0.05  & domain_cevalue < 0.05 & domain_ievalue < 0.05)
pdam.s=pdam.f %>% group_by(query_name) %>% tally()

pver=read_domtblout("Pver_domains")
pver.f=subset(pver, sequence_evalue < 0.05  & domain_cevalue < 0.05 & domain_ievalue < 0.05)
pver.s=pver.f %>% group_by(query_name) %>% tally()

spis=read_domtblout("Spis_domains")
spis.f=subset(spis, sequence_evalue < 0.05  & domain_cevalue < 0.05 & domain_ievalue < 0.05)
spis.s=spis.f %>% group_by(query_name) %>% tally()

# merge results
all=merge(pver.s, pdam.s, by = "query_name", all=T )
all2=merge(all, spis.s, by = "query_name", all=T )
rownames(all2)=all2[,1] 
all2=all2[,-1]
colnames(all2)=c("Pver", "Pdam", "Spis")
all2[is.na(all2)] <- 0
#write.table(all2, "domain_counts.txt", quote = F)

# normalize by total domains found
all.n=sweep(all2,2,colSums(all2),"/")
#write.table(all.n, "domain_counts_norm.txt", quote = F, sep = "\t")

#####################
# Fisher exact test #
#####################

dom.cnt=all2
nor_in = data.frame(Pver = sum(dom.cnt$Pver)-dom.cnt$Pver, Pdam = sum(dom.cnt$Pdam)-dom.cnt$Pdam, Spis = sum(dom.cnt$Spis)-dom.cnt$Spis)

for(i in 1:nrow(dom.cnt)){Pdam.m = matrix(c(dom.cnt$Pver[i],nor_in$Pver[i],dom.cnt$Pdam[i],nor_in$Pdam[i]), nrow=2)}
Pdam.m=lapply(1:nrow(dom.cnt), function(x) Pdam.m = matrix(c(dom.cnt$Pver[x],nor_in$Pver[x],dom.cnt$Pdam[x],nor_in$Pdam[x]), nrow=2))
Pdam.f=lapply(Pdam.m, function(x) fisher.test(x))
dom.cnt$Pdam.pval=lapply(1:nrow(dom.cnt), function(x) Pdam.f[[x]][["p.value"]])
dom.cnt$Pdam.qval=p.adjust(dom.cnt$Pdam.pval, method = "fdr", n=nrow(dom.cnt))


for(i in 1:nrow(dom.cnt)){Spis.m = matrix(c(dom.cnt$Pver[i],nor_in$Pver[i],dom.cnt$Spis[i],nor_in$Spis[i]), nrow=2)}
Spis.m=lapply(1:nrow(dom.cnt), function(x) Spis.m = matrix(c(dom.cnt$Pver[x],nor_in$Pver[x],dom.cnt$Spis[x],nor_in$Spis[x]), nrow=2))
Spis.f=lapply(Spis.m, function(x) fisher.test(x))
dom.cnt$Spis.pval=lapply(1:nrow(dom.cnt), function(x) Spis.f[[x]][["p.value"]])
dom.cnt$Spis.qval=p.adjust(dom.cnt$Spis.pval, method = "fdr", n=nrow(dom.cnt))


message("Number of enriched domains in Pver-Pdam: ", nrow(subset(dom.cnt, Pdam.qval < 0.05))) ## 86
message("Number of enriched domains in Pver-Spis: ", nrow(subset(dom.cnt, Spis.qval < 0.05))) ## 70

enriched_domains=rownames(subset(dom.cnt, Pdam.qval < 0.05 | Spis.qval < 0.05))

###########
# Heatmap #
###########
library(pheatmap)
all.f=subset(all.n, rownames(all.n) %in% enriched_domains)
hp_in=as.matrix(all.f[order(all.f$Pver,decreasing = TRUE),])

pdf("Domain_heatmaps_enriched.pdf", width = 8, height =15, pointsize = 6) 
#heatmap.2(hp_in, scale = "row", key = F, Rowv=F, dendrogram = "none",  density.info = "none", trace = "none",  cexRow=0.5, cexCol=0.8, margins=c(7,7), offsetRow=1, col=grey(seq(1,0,-0.01)))
#heatmap.2(hp_in, scale = "row", key = F, Rowv=T, dendrogram = "none",  density.info = "none", trace = "none",  cexRow=0.5, cexCol=0.8, margins=c(7,7), offsetRow=1, col=colorRampPalette(c("navy", "white", "firebrick3"))(50))
out <-pheatmap(hp_in, color = colorRampPalette(c("navy", "white", "firebrick3"))(50), angle_col = "0",  scale = "row", cluster_col = FALSE)
dev.off()
