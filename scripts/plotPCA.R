#!/usr/bin/Rscript
# Usage: Rscript -i infile.covar -c component1-component2 -a annotation.file -o outfile.eps


library(optparse)
library(ggplot2)

option_list <- list(make_option(c('-i','--in_file'), action='store', type='character', default=NULL, help='Input file (output from ngsCovar)'),
                    make_option(c('-c','--comp'), action='store', type='character', default=1-2, help='Components to plot'),
                    make_option(c('-a','--annot_file'), action='store', type='character', default=NULL, help='Annotation file with individual classification (2 column TSV with ID and ANNOTATION)'),
                    make_option(c('-o','--out_file'), action='store', type='character', default=NULL, help='Output file')
                    )
opt <- parse_args(OptionParser(option_list = option_list))

# Annotation file is in plink cluster format

#################################################################################

# Read input file
covar <- read.table(opt$in_file, stringsAsFact=F);

# Read annot file
annot <- read.table(opt$annot_file, sep=" ", header=T); # note that plink cluster files are usually tab-separated instead

# Parse components to analyze
comp <- as.numeric(strsplit(opt$comp, "-", fixed=TRUE)[[1]])



# Eigenvalues
eig <- eigen(covar, symm=TRUE);
eig$val <- eig$val/sum(eig$val);
cat(signif(eig$val, digits=3)*100,"\n");


# Write eigenvalues
#write.table(eig, file = "eigen_scores.txt", quote = FALSE)


# Plot
PC <- as.data.frame(eig$vectors)
colnames(PC) <- gsub("V", "PC", colnames(PC))
PC$Color <- factor(annot$CLUSTER)
PC$Alt <- factor(annot$IID)
PC$ID <- factor(annot$FID)


# Write PC components
write.table(PC, file = "PC_scores.txt", quote = FALSE)


#title <- paste("PC",comp[1]," (",signif(eig$val[comp[1]], digits=3)*100,"%)"," / PC",comp[2]," (",signif(eig$val[comp[2]], digits=3)*100,"%)",sep="",collapse="")

#x_axis = paste("PC",comp[1],sep="")
#y_axis = paste("PC",comp[2],sep="")

isoPalette=c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise", "green1", "yellow4", "yellow3", "darkorange4", "brown") ### 25 Colors

### Basic plot ###
#ggplot() + geom_point(data=PC, aes_string(x=x_axis, y=y_axis, color="Pop")) + ggtitle(title) + scale_colour_manual(values=isoPalette)

### Basic plot with labeled individuals ###
#ggplot() + geom_point(data=PC, aes_string(x=x_axis, y=y_axis, color="Pop")) + ggtitle(title) + scale_colour_manual(values=isoPalette) + geom_text(data=PC, aes_string(x=x_axis, y=y_axis, label="Lab", vjust= -0.784), cex=2.5, hjust=-0.3)

### Multi group plot with labeled individuals ###
#PC$Rate <- factor(PC$Rate, levels = c("Low", "Mid", "High"))
#ggplot() + geom_point(data=PC, aes_string(x=x_axis, y=y_axis, color="Rate", shape="Pop")) + ggtitle(title) + scale_colour_manual(values=c("#e6b400", "#e47200", "red")) + scale_shape_manual(values = c(15,17,19,23,14,25,8,9,11,12,7))

# Calculate axis labels with percentages
x_axis_label <- paste0("PC", comp[1], " (", signif(eig$val[comp[1]], digits = 3) * 100, "%)")
y_axis_label <- paste0("PC", comp[2], " (", signif(eig$val[comp[2]], digits = 3) * 100, "%)")

# Create the plot
#ggplot() +
#  geom_point(data = PC, aes_string(x = paste0("PC", comp[1]), y = paste0("PC", comp[2]), color = "Color", shape = "Alt"), size = 5) +  # Increase size of shapes
#  scale_colour_manual(values = c("#000000", "#009E73")) +  # Black and green (colorblind-friendly)
#  scale_shape_manual(values = c(15, 16, 17, 18, 0, 2, 3, 4, 5, 6, 8)) +  # First 4 solid, last 7 unfilled and distinct
#  labs(x = x_axis_label, y = y_axis_label) +  # Set custom axis labels
#  theme_classic()  # Optional, for a cleaner plot background

ggplot() +
  geom_point(data = PC, aes_string(x = paste0("PC", comp[1]), y = paste0("PC", comp[2]), color = "Color", shape = "Alt"), size = 5) +  # Adjust point size for readability
  scale_colour_manual(values = c("#000000", "#009E73")) +  # Black and green (colorblind-friendly)
  scale_shape_manual(values = c(15, 16, 17, 18, 0, 2, 3, 4, 5, 6, 8)) +  # First 4 solid, last 7 unfilled and distinct
  labs(x = x_axis_label, y = y_axis_label) +  # Set custom axis labels
  theme_classic(base_size = 10) +  # Set base text size for improved readability
  theme(
    axis.text = element_text(size = 10),       # Axis tick labels size
    axis.title = element_text(size = 10),    # Axis titles size
    legend.text = element_text(size = 10),    # Legend text size
    legend.title = element_text(size = 10),   # Legend title size
    plot.margin = margin(4, 4, 4, 4, "pt")   # Adjust margins to avoid cutoff
  ) +
  coord_fixed(ratio = 1.5)  # Adjust aspect ratio to make the plot rectangular

#ggsave(opt$out_file, width = 5.5, height = 8.25, units = "cm", dpi = 300)
ggsave(opt$out_file)
unlink("Rplots.pdf", force=TRUE)
