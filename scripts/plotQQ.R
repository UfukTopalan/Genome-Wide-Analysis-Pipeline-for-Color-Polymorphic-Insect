setwd("/home/ismail/Research/isophya/PCA/NEW/")

# read loadings
loadings <- read.table("isophya71.pca.loadings", header = TRUE)

# Function to draw the figure (reused for PDF and PNG)
draw_qq_plots <- function() {
  
  par(
    mfrow = c(2, 1),
    mar   = c(5, 5, 2.5, 2),
    cex.axis = 1.1,
    cex.lab  = 1.2
  )
  
  # --- Panel A: PC1 ---
  qqnorm(
    loadings$PC1,
    pch  = 19,
    cex  = 0.6,
    col  = adjustcolor("#00441B", alpha.f = 0.6),
    xlab = "Theoretical quantiles",
    ylab = "PC1 loadings"
  )
  qqline(loadings$PC1, lwd = 2)
  mtext("A", side = 3, line = 0.5, adj = 0, cex = 1.3, font = 2)
  
  # --- Panel B: PC2 ---
  qqnorm(
    loadings$PC2,
    pch  = 19,
    cex  = 0.6,
    col  = adjustcolor("#00441B", alpha.f = 0.6),
    xlab = "Theoretical quantiles",
    ylab = "PC2 loadings"
  )
  qqline(loadings$PC2, lwd = 2)
  mtext("B", side = 3, line = 0.5, adj = 0, cex = 1.3, font = 2)
}

# -------------------
# Save as PDF (A4)
# -------------------
pdf(
  file = "QQplots_PC1_PC2.pdf",
  width = 8.27,   # A4 width (inches)
  height = 11.69  # A4 height (inches)
)
draw_qq_plots()
dev.off()

# -------------------
# Save as PNG (A4, high resolution)
# -------------------
png(
  filename = "QQplots_PC1_PC2.png",
  width  = 2480,  # A4 width at 300 dpi
  height = 3508,  # A4 height at 300 dpi
  res    = 300
)
draw_qq_plots()
dev.off()
