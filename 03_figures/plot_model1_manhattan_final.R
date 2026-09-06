args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]
output_prefix <- args[2]

dat <- read.table(input_file, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, check.names = FALSE)

dat$chr <- as.numeric(dat$chr)
dat$CHR_POS <- as.numeric(dat$CHR_POS)
dat$UNADJ <- as.numeric(dat$UNADJ)

dat <- dat[complete.cases(dat[, c("chr", "CHR_POS", "UNADJ")]), ]
dat <- dat[dat$UNADJ > 0, ]

dat$logP <- -log10(dat$UNADJ)

chr_max <- tapply(dat$CHR_POS, dat$chr, max)
offset <- c(0, cumsum(chr_max)[-length(chr_max)])
names(offset) <- names(chr_max)

dat$x <- dat$CHR_POS + offset[as.character(dat$chr)]

chr_mid <- offset + chr_max / 2

bonf <- 0.05 / 1430007

plot_manhattan <- function(outfile, type = "pdf") {

  if (type == "pdf") {
    pdf(outfile, width = 11, height = 6)
  } else {
    png(outfile, width = 3300, height = 1800, res = 300)
  }

  plot(dat$x, dat$logP,
       pch = 16,
       cex = 0.35,
       col = ifelse(dat$chr %% 2 == 0, "grey50", "grey20"),
       xaxt = "n",
       xlab = "Daphnia magna chromosome",
       ylab = expression(-log[10](P)),
       main = "Model 1 co-GWAS")

  axis(1, at = chr_mid, labels = names(chr_max), tick = FALSE)

  abline(h = -log10(bonf),
         col = "red",
         lty = 2,
         lwd = 1.5)

  dev.off()
}

plot_manhattan(paste0(output_prefix, ".pdf"), "pdf")
plot_manhattan(paste0(output_prefix, ".png"), "png")
