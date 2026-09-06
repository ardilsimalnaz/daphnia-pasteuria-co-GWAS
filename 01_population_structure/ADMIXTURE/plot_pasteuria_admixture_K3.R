 q <- read.table("pasteuria.admix.3.Q")
fam <- read.table("pasteuria.admix.fam")
groups <- read.table("pasteuria_lineages.txt", header=TRUE, sep="\t")

group <- groups$group[match(fam$V2, groups$ID)]

order_groups <- c("Alpha-1","Alpha-2","Alpha-3","Beta","Gamma","Mixed")
ord <- order(factor(group, levels=order_groups))

q <- q[ord,]
group <- group[ord]

counts <- table(factor(group, levels=order_groups))
ends <- cumsum(counts)

png("Pasteuria_ADMIXTURE_K3_final.png", width=1600, height=700, res=150)

par(mar=c(5,5,2,9), xpd=NA)

bp <- barplot(t(as.matrix(q)),
              col=c("red","gold","blue"),
              border=NA,
              space=0,
              axes=FALSE,
names.arg=rep("", nrow(q)),
              ylab="Ancestry proportion")

axis(2)

mids <- sapply(1:6, function(i) mean(bp[group == order_groups[i]]))

axis(1, at=mids, labels=order_groups,
     tick=FALSE, cex.axis=0.8)

for(i in ends[-6]) {
    segments((bp[i]+bp[i+1])/2, 0,
             (bp[i]+bp[i+1])/2, 1,
             lty=2)
}

legend("topright",
       inset=c(-0.18,0),
       legend=paste0("Component ",1:3),
       fill=c("red","gold","blue"),
       bty="n")

dev.off()
