# Genome size estimation using Jellyfish

setwd("~/Genomes/corals/Pocillopora.verrucosa/05_Pver_lib4_genome_properties")

Pver_31mer <- read.table("Pver_31mer.histo") #25622 lanes
# het_peak 76 4777149
# hom_peak 154 3760129

# (1)
sum(as.numeric(Pver_31mer[23:116,1]*Pver_31mer[23:116,2]))/76
# 160,624,609

# (2)
sum(as.numeric(Pver_31mer[117:25622,1]*Pver_31mer[117:25622,2]))/154
# 260,217,036

# (3) summing those values.
160624609 + 260217036 

# 420,841,645 Diploid genome size 31mer

Pver_25mer <- read.table("Pver_25mer.histo") #29660 lanes
# het_peak 79 3985807
# hom_peak 159 3676856

# (1)
sum(as.numeric(Pver_25mer[23:116,1]*Pver_25mer[23:116,2]))/79
# 134,341,370

# (2)
sum(as.numeric(Pver_25mer[117:29660,1]*Pver_25mer[117:29660,2]))/159
# 272,892,790

# (3) summing those values.
134341370 + 272892790
# 407,234,160  Diploid genome size size 25mer

marks <- c(0,1000000,2000000,3000000,4000000,5000000)
marks_label <- c(0,1,2,3,4,5)

png("Khmer_spectra_jellyfish.png", width = 600, height = 300) 
# Create a plot
par(mar = c(5, 7, 2, 2))
plot(Pver_31mer[4:300,], type="l", col="gray18", lwd = 2, lty=4, axes = F, frame.plot=F, ann=F)
axis(1, pos=0, las=1)
title(xlab = "Kmer-Coverage", cex.lab = 1, line = 2)
axis(2, pos=0, las=1, at=marks, labels = marks_label)
title(ylab = "Kmer counts frequency (millions)", cex.lab = 1, line = 1.5)
lines(Pver_25mer[4:300,], type="l", col="gray18", las=1, lwd = 2)
legend("topright", c("k = 25", "k = 31"), cex=1, col = "gray18", lty=c(1,4) , bty = "n")
dev.off() 
