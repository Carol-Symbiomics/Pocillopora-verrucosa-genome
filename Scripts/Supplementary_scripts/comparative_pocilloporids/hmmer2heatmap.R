
library(rhmmer)
library(gplots)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/Pver_genome/")

# import hmmer results
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
write.table(all2, "domain_counts.txt", quote = F)

# plot heatmap (most abundant 100 domains)
all_low_Var=all2[apply(all2, 1, var) != 0,] #remove rows with no variation
topDomains=rownames(all_low_Var[order(rowSums(all_low_Var),decreasing = TRUE),][1:100,])
hp_in=as.matrix(subset(all2, rownames(all2) %in% topDomains))

pdf("Domain_heatmaps.pdf",  width = 3, height =10, pointsize = 12) 
heatmap.2(hp_in, scale = "row", key = TRUE, Colv=FALSE, dendrogram = "row",  density.info = "none", trace = "none", cexRow=0.3, cexCol=0.8, margins=c(7,7), offsetRow=1, col=grey(seq(1,0,-0.01)))
dev.off()


