q <- read.table("daphnia_faststructure.6.meanQ")

png("Daphnia_fastStructure_K6.png",
    width=1600, height=700, res=150)

par(mar=c(5,5,2,9), xpd=NA)

barplot(t(as.matrix(q)),
        col=rainbow(6),
        border=NA,
        space=0,
        names.arg=rep("", nrow(q)),
        xlab="Daphnia samples",
        ylab="Ancestry proportion")

legend("topright",
       inset=c(-0.18,0),
       legend=paste0("Component ",1:6),
       fill=rainbow(6),
       bty="n")

dev.off()
