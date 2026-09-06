pca <- read.table("daphnia.pca.eigenvec",
                  header=TRUE, sep="\t", comment.char="")

png("Daphnia_PCA_paperstyle.png",
    width=1000, height=1000, res=150)

plot(pca$PC1, pca$PC2,
     pch=21,
     bg=rgb(0.2, 0.6, 1, 0.55),
     col="grey30",
     cex=1.3,
     xlab="PC 1 (2.86% of variance)",
     ylab="PC 2 (2.16% of variance)",
     bty="o",
     las=1)

dev.off()
