# Run this as an R script

library(showtext)
library(dplyr)
library(ggplot2)
library(ape)
showtext_auto()
library(viridis)
library(scales)

workdir <- "/mnt/storage11/sophie/david/filtered_vcfs/pca" # Working directory with plink files
prefix <- "pfal_wgs_and_amp" # Prefix for plink files
metadata <- "/mnt/storage11/sophie/david/wgs_amp_metadata.csv" # File path to metadata

calc_variance_explained <- function(pc_points) {
    vars <- round(pc_points$eig / sum(pc_points$eig) * 100, 1)
    names(vars) <- paste0("PC", seq_len(length(vars)))
    vars
}

# METADATA
met <- read.table(metadata, sep = ",", stringsAsFactors = FALSE, header = TRUE)

#### DIST#
dist <- read.table(file.path(workdir, paste0(prefix, ".dist")), header = FALSE)
id <- read.table(file.path(workdir, paste0(prefix, ".dist.id")))

desc <- id %>% left_join(met, by = c("V1" = "sample"))

dist_m <- as.matrix(dist)
colnames(dist_m) <- desc$V1
rownames(dist_m) <- desc$V1

### Debugging NA distances ###
na_pos <- which(is.na(dist_m), arr.ind = TRUE)
if (nrow(na_pos) == 0) {
  cat("No NA values in the distance matrix.\n")
} else {
  cat("Found NA distances between sample pairs:\n")
  apply(na_pos, 1, function(idx) {
    cat(sprintf("  %s <-> %s\n", rownames(dist_m)[idx[1]], colnames(dist_m)[idx[2]]))
  })
  bad_samples <- unique(c(rownames(dist_m)[na_pos[,1]], colnames(dist_m)[na_pos[,2]]))
  cat("\nSamples with any NA distances:\n")
  cat(sprintf("  %s\n", bad_samples), sep = "")
}

# Impute NA distances to the median
med <- median(dist_m, na.rm = TRUE)
med <- 0
dist_m[is.na(dist_m)] <- med

# PCA (MDS)
cmd <- cmdscale(dist_m, k = 10, eig = TRUE, x.ret = TRUE)
vars <- calc_variance_explained(cmd)

# Prepare data frame for plotting
df <- as.data.frame(cmd$points, stringsAsFactors = FALSE)
df$sample <- rownames(df)
df$sequencing_method <- gsub("_", " ", desc$sequencing_method)
df$day <- gsub("_", " ", desc$day)

# Reorder so 'sample' is first
df <- df[, c("sample", setdiff(names(df), "sample"))]
# Rename PCA columns
colnames(df)[-c(1, ncol(df)-1, ncol(df))] <- gsub("^V", "PC", colnames(df)[-c(1, ncol(df)-1, ncol(df))])

color_by <- "sequencing_method"
my_colours <- c("WGS" = "#7678ed", "Amplicon" = "darkorange2")

# Plot without labels, include legend
png("WGS_vs_Amp_data.png", width = 800, height = 600)
plt <- ggplot(data = df, aes(x = PC1, y = PC2, color = !!sym(color_by))) +
    geom_point(size = 3) +
    labs(
      x = paste0("PC1 (", vars["PC1"], "%)"),
      y = paste0("PC2 (", vars["PC2"], "%)"),
      title = "Amplicon vs WGS samples",
      color = "Sequencing Method"
    ) +
    scale_color_manual(values = my_colours) +
    scale_x_continuous(labels = label_number()) +
    scale_y_continuous(labels = label_number()) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      plot.title = element_text(hjust = 0.5),
      plot.margin = margin(t = 10, r = 40, b = 30, l = 10, unit = "pt")
    )
print(plt)
dev.off()
