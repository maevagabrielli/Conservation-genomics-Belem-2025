### Barplot ROH length

# =========================
# Setup
# =========================
setwd("/Users/mgabrielli/Documents/Postdoc/Data/")
meta <- read.table("Podarcis_inds_sp_pop.txt", header = FALSE)
colnames(meta) <- c("ind", "species", "pop")

# Individuals to plot
meta <- meta[meta$ind %in% c("LC03", "ST02", "DB23461"), ]
genome_size <- 42074563  # genome length (bp), here of chromosome 17 (total genome size: 1.43Mb)

# ROH length classes
roh_classes <- list(c(5e5, 1e6), c(1e6, 2e6), c(2e6, 5e6), c(5e6, Inf))
class_names <- c("500Kb-1Mb (100 gen)","1Mb-2Mb (50 gen)","2Mb-5Mb (25 gen)",">5Mb (10 gen)")

# =========================
# Compute ROH proportions
# =========================
tab <- matrix(0, nrow = 4, ncol = nrow(meta))
colnames(tab) <- meta$ind
rownames(tab) <- class_names

setwd("/Users/mgabrielli/Documents/Postdoc/Projects/Brazil/Belem2025/workshop/")
for (j in seq_len(nrow(meta))) {
  file <- paste0("length_ROH_", meta$ind[j], ".txt")
  if (!file.exists(file) || file.info(file)$size == 0) next # check that file exists and is not empty (if no ROH)
  roh <- read.table(file)[,1]
  for (k in 1:4) { # loop over the 4 ROH classes
    tab[k, j] <- sum(roh[roh >= roh_classes[[k]][1] & roh < roh_classes[[k]][2]]) / genome_size
  }
}

# =========================
# Plot and save in a PDF file
# =========================

pdf("Barplot_length_ROH.pdf", width = 12, height = 7)

barplot(tab,
  col = colors()[c(24,12,89,23)],
  border = "white",
  space = 0.05,
  las = 2,
  ylim = c(0, 1),
  ylab = "Proportion of genome in ROH",
  main = "ROH length distribution per individual",
  legend.text = rownames(tab))

dev.off()