### This whole thing is the Lhat function
### Everything happens inside of here. 
myMonocle <- function(sce, metadata) {
  
  # Make a monocle object
  X <- t(counts(sce))
  gene_meta <- data.frame(gene_short_name = rownames(t(X)))
  rownames(gene_meta) <- rownames(t(X))
  cds <- new_cell_data_set(t(X),
                           cell_metadata = metadata, 
                           gene_metadata = gene_meta
  )
  
  # Compute the high variance genes.
  sizeFactors(sce) <- librarySizeFactors(sce)
  sce <- logNormCounts(sce)
  dec.cm<- modelGeneVar(sce)
  hvg.cm.var <- getTopHVGs(dec.cm, fdr.threshold=0.005)
  if (length(hvg.cm.var) > 2500) {
    hvg.cm.var <- hvg.cm.var[1:2500] 
  }
  
  ### I believe this does size factor norm, log(X+1), and subsets to HVGs. 
  cds <- preprocess_cds(cds, num_dim = 100, use_genes = hvg.cm.var)
  
  sce <- estimate_cycle_position(sce, species="human", gname.type="SYMBOL")
  colData(cds)$phase <- colData(sce)$tricyclePosition
  
  ### I think that this combination of residual model formula works well and seems reasonable. 
  ### The one thing is that Scores should be calculated for THIS set only. 
  cds <- align_cds(cds,alignment_group = "individual", 
                   residual_model_formula_str = "~ phase + percent.mito")
  cds <- reduce_dimension(cds)
  
  ## You are required to run this before running pseudotime.
  cds <- cluster_cells(cds, partition_qval = 0.05)
  
  ## I have use_partition=FALSE because I want one single trajectory!
  ## I don't want multiple roots. 
  cds <- learn_graph(cds, use_partition=FALSE)
  
  #### This is how I pick the IPSC cell as a root.
  #### Copied from monocle documentation. 
  cell_ids <- which(colData(cds)[, "type"] == "IPSC")
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  ### Finally!! Get the pseudotime. 
  cds <- order_cells(cds, root_pr_nodes=root_pr_nodes)
  return(cds)
}