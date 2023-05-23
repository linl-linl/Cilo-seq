library(grid)
library(futile.logger)
library(VennDiagram)
library(grDevices)
bulk <- read.table("bulk_2000000_genelist_FPKM.txt", header = T)
seq <- read.table("seq_2000000_genelist_FPKM.txt", header = T)
seq2 <- read.table("seq2_2000000_genelist_FPKM.txt", header = T)
tube<- read.table("tube_2000000_genelist_FPKM.txt", header = T)
venn.plot1 <- venn.diagram(list("Digital-RNA-Seq"=seq$GENE, "Bulk"=bulk$GENE), filename=NULL, 
                           height = 5000, width = 5000, lwd = 3, lty = 1, 
                           fill = c("#038FF4","#A0A0A0"), cex = 2, fontface = "bold", 
                           cat.cex = 2, cat.fontface = "bold", cat.pos = c(25, -35), 
                           cat.dist = c(0.055, 0.045), rotation.degree = 180)
venn.plot2 <- venn.diagram(list("Digital-RNA-Seq2"=seq2$GENE, "Bulk"=bulk$GENE), filename=NULL, 
                           height = 5000, width = 5000, lwd = 3, lty = 1, 
                           fill = c("#FF390C","#A0A0A0"), cex = 2, fontface = "bold", 
                           cat.cex = 2, cat.fontface = "bold", cat.pos = c(24, -35), 
                           cat.dist = c(0.05, 0.045), rotation.degree = 180)
venn.plot3 <- venn.diagram(list("Tube"=tube$GENE, "Bulk"=bulk$GENE), filename=NULL, 
                           height = 5000, width = 5000, lwd = 3, lty = 1, 
                           fill = c("#FFB402","#A0A0A0"), cex = 2, fontface = "bold", 
                           cat.cex = 2, cat.fontface = "bold", cat.pos = c(24, -35), 
                           cat.dist = c(0.045, 0.045), rotation.degree = 180)
venn.plot4 <- venn.diagram(list("Digital-RNA-Seq"=seq$GENE, "Digital-RNA-Seq2"=seq2$GENE), filename=NULL, 
                           height = 5000, width = 5000, lwd = 3, lty = 1, 
                           fill = c("#038FF4","#FF390C"), cex = 2, fontface = "bold", 
                           cat.cex = 2, cat.fontface = "bold", cat.pos = c(25, -20), 
                           cat.dist = c(0.055, 0.045), rotation.degree = 180)
venn.plot5 <- venn.diagram(list("Tube"=tube$GENE, "Digital-RNA-Seq2"=seq2$GENE), filename=NULL, 
                           height = 5000, width = 5000, lwd = 3, lty = 1, 
                           fill = c("#FFB402","#FF390C"), cex = 2, fontface = "bold", 
                           cat.cex = 2, cat.fontface = "bold", cat.pos = c(25, -20), 
                           cat.dist = c(0.055, 0.045), rotation.degree = 180)
venn.plot6 <- venn.diagram(list("Digital-RNA-Seq"=seq$GENE, "Tube"=tube$GENE), filename=NULL, 
                           height = 5000, width = 5000, lwd = 3, lty = 1, 
                           fill = c("#038FF4","#FFB402"), cex = 2, fontface = "bold", 
                           cat.cex = 2, cat.fontface = "bold", cat.pos = c(25, -20), 
                           cat.dist = c(0.055, 0.045), rotation.degree = 180)
pdf(file="seq_bulk.pdf")
grid.draw(venn.plot1)
dev.off()
pdf(file="seq2_bulk.pdf")
grid.draw(venn.plot2)
dev.off()
pdf(file="tube_bulk.pdf")
grid.draw(venn.plot3)
dev.off()
pdf(file="seq_seq2.pdf")
grid.draw(venn.plot4)
dev.off()
pdf(file="tube_seq2.pdf")
grid.draw(venn.plot5)
dev.off()
pdf(file="seq_tube.pdf")
grid.draw(venn.plot6)
dev.off()
