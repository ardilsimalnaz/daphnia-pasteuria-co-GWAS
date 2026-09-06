pca <- read.table("pasteuria.pca.all.eigenvec",
                  header=TRUE, sep="\t", comment.char="")
eig <- scan("pasteuria.pca.all.eigenval")
groups <- read.table("pasteuria_lineages.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

pc1var <- round(100 * eig[1] / sum(eig), 2)
pc2var <- round(100 * eig[2] / sum(eig), 2)

dat <- merge(pca, groups, by.x="IID", by.y="ID")

png("Pasteuria_PCA_subgroups.png", width=1200, height=1000, res=150)
par(mar=c(5, 7, 2, 2))

plot(dat$PC1, dat$PC2,
     type="n",
     xlim=c(min(dat$PC1)-0.01, max(dat$PC1)+0.01),
     ylim=c(min(dat$PC2)-0.01, max(dat$PC2)+0.03),
     xlab=paste0("PC 1 (", pc1var, "% of variance)"),
     ylab="",
     cex.lab=1.4,
     cex.axis=1.3,
     las=1)
mtext(paste0("PC 2 (", pc2var, "% of variance)"),
      side=2, line=4.5, cex=1.4)

cols <- c(
  "Alpha-1" = "tomato",
  "Alpha-2" = "orange",
  "Alpha-3" = "gold",
  "Beta"    = "limegreen",
  "Gamma"   = "dodgerblue",
  "Mixed"   = "white"
)

for(g in c("Mixed","Alpha-1","Alpha-2","Alpha-3","Beta","Gamma")) {
  sub <- dat[dat$group == g, ]
  points(sub$PC1, sub$PC2,
         pch=21,
         bg=cols[g],
         col="grey30",
         cex=1.4)
}

alpha <- dat$group %in% c("Alpha-1","Alpha-2","Alpha-3")
beta  <- dat$group == "Beta"
gamma <- dat$group == "Gamma"

text(mean(dat$PC1[alpha]), mean(dat$PC2[alpha]) + 0.015, "Alpha", cex=1.2)
text(mean(dat$PC1[beta]) + 0.01, mean(dat$PC2[beta]) + 0.015, "Beta", cex=1.2)
text(mean(dat$PC1[gamma]) + 0.01, mean(dat$PC2[gamma]) + 0.015, "Gamma", cex=1.2)

legend("topright",
       legend=c("Alpha-1","Alpha-2","Alpha-3","Beta","Gamma","Mixed"),
       pt.bg=c("tomato","orange","gold","limegreen","dodgerblue","white"),
       pch=21,
       pt.cex=1.2,
       col="grey30",
       bty="n",
       cex=1.1)

dev.off()
