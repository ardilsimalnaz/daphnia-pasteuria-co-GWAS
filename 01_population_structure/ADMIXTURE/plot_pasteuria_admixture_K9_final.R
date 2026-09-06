q <- read.table("pasteuria.admix.9.Q")
fam <- read.table("pasteuria.admix.fam",
                  stringsAsFactors=FALSE)

groups <- read.table("pasteuria_lineages.txt",
                     header=TRUE,
                     sep="\t",
                     stringsAsFactors=FALSE)

group <- groups$group[match(fam$V2, groups$ID)]
stopifnot(nrow(q) == nrow(fam))
stopifnot(!any(is.na(group)))

levels_order <- c("Alpha-1", "Alpha-2", "Alpha-3",
                  "Beta", "Gamma", "Mixed")

ord <- order(factor(group, levels=levels_order))

q <- q[ord, , drop=FALSE]
group <- group[ord]
counts <- sapply(levels_order,
                 function(g) sum(group == g))

starts <- c(1, head(cumsum(counts), -1) + 1)
ends <- cumsum(counts)

png("Pasteuria_ADMIXTURE_K9_FINAL.png",
    width=1700, height=750, res=150)

par(mar=c(7,5,2,11), xpd=TRUE)

bp <- barplot(t(as.matrix(q)),
              col=rainbow(9),
              border=NA,
              space=0,
              xaxt="n",
              xlab="",
              ylab="Ancestry proportion")

mids <- (bp[starts] + bp[ends]) / 2

axis(1,
     at=mids,
     labels=levels_order,
     tick=FALSE,
     line=1,
     cex.axis=0.85,
gap.axis=-1)
for(i in ends[-length(ends)]) {
    abline(v=(bp[i] + bp[i+1])/2,
           lty=2)
}

legend("topright",
       inset=c(-0.19,0),
       legend=paste0("Component ", 1:9),
       fill=rainbow(9),
       bty="n",
       cex=0.75)

dev.off()

print(counts)
