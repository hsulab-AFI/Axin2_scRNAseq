# Trajectory_analysis_using Monocle3
# load packages
library(Seurat)
library(SingleCellExperiment)
remotes::install_github("satijalab/seurat-wrappers")
library(SeuratWrappers)
library(monocle3)
library(ggplot2)
library(dplyr)
library(Matrix)
library(tidyverse)
library(cowplot)
library(patchwork)

# Set up monocle3 cell_data_set object using the SeuratWrappers
cds <- as.cell_data_set(Axin2.combined)
cds <- cluster_cells(cds, resolution=1e-6)

p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

# Subsetting partitions
integrated.sub <- subset(as.Seurat(cds, assay = NULL), monocle3_partitions == 1)
cds <- as.cell_data_set(integrated.sub)

# Trajectory analysis
cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE)

plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size = 1.5)

# Color cells by pseudotime
#cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 2]))
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           label_roots = TRUE,
           trajectory_graph_color = "grey60")&
  theme(axis.text=element_text(size=18),axis.title=element_text(size=18))

integrated.sub <- as.Seurat(cds, assay = NULL)
FeaturePlot(integrated.sub, "monocle3_pseudotime")&
  theme(axis.text=element_text(size=18),axis.title=element_text(size=18))


# Identify genes that change as a function of pseudotime
cds_graph_test_results <- graph_test(cds,
                                     neighbor_graph = "principal_graph",
                                     cores = 8)
rowData(cds)$gene_short_name <- row.names(rowData(cds))

head(cds_graph_test_results, error=FALSE, message=FALSE, warning=FALSE)

deg_ids <- rownames(subset(cds_graph_test_results[order(cds_graph_test_results$morans_I, decreasing = TRUE),], q_value < 0.05))

plot_cells(cds,
           genes=head(deg_ids),
           show_trajectory_graph = FALSE,
           label_cell_groups = FALSE,
           label_leaves = FALSE)
