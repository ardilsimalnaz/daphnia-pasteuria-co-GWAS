args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]
output_prefix <- args[2]

suppressPackageStartupMessages(library(circlize))

dat <- read.delim(input_file)
map <- read.csv("contig2chrom.csv")

dat$LOG10_P <- -log10(dat$UNADJ)

dat$daphBase <- dat$DAPH_POS / 125578570
dat$pastBase <- dat$PAST_POS / 1723598

chr_start <- tapply(map$full_position, map$chr, min)
chr_end <- tapply(map$full_position + map$contig_length, map$chr, max)
chr_mid <- (chr_start + chr_end) / 2 / 125578570

n <- nrow(dat)

genome <- factor(
  c(rep("Host", n), rep("Parasite", n)),
  levels = c("Host", "Parasite")
)

x <- c(dat$daphBase, dat$pastBase)
y <- c(dat$LOG10_P, dat$LOG10_P)

plotdata <- data.frame(genome, x, y)

plot_circle <- function(outfile, type) {

  if (type == "pdf") {
    pdf(outfile, width = 7, height = 7)
  } else {
    png(outfile, width = 2100, height = 2100, res = 300)
  }

  par(mar = c(2, 2, 4, 2))

  circos.clear()

  circos.par(
    gap.degree = 15,
    start.degree = -8,
    clock.wise = TRUE,
    cell.padding = c(0, 0, 0, 0)
  )

  circos.initialize(
    factors = c("Host", "Parasite"),
    xlim = cbind(c(0, 0), c(1, 1))
  )

  circos.trackPlotRegion(
    factors = plotdata$genome,
    x = plotdata$x,
    y = plotdata$y,
    track.height = 0.16,
    bg.col = "grey90",
    bg.border = "grey40",
    panel.fun = function(x, y) {

      sector <- CELL_META$sector.index

      circos.points(x, y, pch = 16, cex = 0.22)

      if (sector == "Parasite") {

        mb <- c(0, 0.4, 0.8, 1.2, 1.6)

        circos.axis(
          h = "top",
          major.at = mb / 1.723598,
          labels = mb,
          labels.cex = 0.55,
          minor.ticks = 1,
          labels.niceFacing = TRUE
        )
      }

      if (sector == "Host") {

        circos.axis(
          h = "top",
          major.at = chr_mid,
          labels = names(chr_mid),
          labels.cex = 0.6,
          minor.ticks = 0,
          labels.niceFacing = TRUE
        )
      }
    }
  )

  link_col <- adjustcolor("firebrick", alpha.f = 0.12)

  for (i in 1:n) {
    circos.link(
      "Host", dat$daphBase[i],
      "Parasite", dat$pastBase[i],
      h = 1,
      col = link_col,
      border = NA,
      lwd = 0.5
    )
  }

  title("Co-GWAS Model 1", cex.main = 1.5)

  text(
    0, 1.12,
    expression(italic(P.~ramosa) ~ "genome"),
    cex = 1.1,
    xpd = NA
  )

  text(
    0, -1.12,
    expression(italic(D.~magna) ~ "genome"),
    cex = 1.1,
    xpd = NA
  )

  dev.off()
}

plot_circle(paste0(output_prefix, ".pdf"), "pdf")
plot_circle(paste0(output_prefix, ".png"), "png")

cat("Associations plotted:", n, "\n")
