q <- read.table("daphnia.admix_num.6.Q")

png("Daphnia_ADMIXTURE_K6_FINAL.png",
    width=1600, height=700, res=150)

par(mar=c(5, 5, 2, 9), xpd=TRUE)

barplot(t(as.matrix(q)),
        col=rainbow(6),
        border=NA,
        space=0,
        xlab="Daphnia samples",
        ylab="Ancestry proportion")

legend("topright",
       inset=c(-0.22, 0),
       legend=paste0("Component ", 1:6),
       fill=rainbow(6),
       bty="n",
       cex=0.8)

dev.off()
