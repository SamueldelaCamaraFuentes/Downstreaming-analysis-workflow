interactions <- function(df){
  #Subset de up-regulated y down-regulated
  
  up_regulated <- subset(df, df$expression == "Up-regulated")
  down_regulated <- subset(df, df$expression == "Down-regulated")
  
  up_mapped <- string_db$map(up_regulated, "Protein", removeUnmappedRows = TRUE)
  
  down_mapped <- string_db$map(down_regulated, "Protein", removeUnmappedRows = TRUE)
  
  par(mfrow=c(1,1))
  hits_up <- up_mapped$STRING_id[1:100] 
  hits_down <- down_mapped$STRING_id[1:100]
  filename = "Interaction network of upregulated proteins.tiff"
  tiff(filename = filename, units="in", width=9, height=9, res=300)
  string_db$plot_network(hits_up)
  dev.off()
  filename = "Interaction network of downregulated proteins.tiff"
  tiff(filename = filename, units="in", width=9, height=9, res=300)
  string_db$plot_network(hits_down)
  dev.off()
  
  interactions <- list(up_mapped, down_mapped)
  return(interactions)
  
}


##################################################################################
#igraph analysis 

igraph_analysis <- function(interactions, taxonid, score) {
  string_db <- STRINGdb$new(version = "11.5", species = taxonid, score_threshold = score, input_directory = "", protocol = "http")
  
  hits_upregulated <- interactions[[1]]
  hits_downregulated <- interactions[[2]]
  
  subgraph_up_proteins <- string_db$get_subnetwork(hits_upregulated$STRING_id[1:100])
  subgraph_down_proteins <- string_db$get_subnetwork(hits_downregulated$STRING_id[1:100])
  
  graph_analysis_up <- calculate_graph_measures(subgraph_up_proteins, hits_upregulated)
  graph_analysis_down <- calculate_graph_measures(subgraph_down_proteins, hits_downregulated)
  
  list(Upregulated = graph_analysis_up, Downregulated = graph_analysis_down)
}

calculate_graph_measures <- function(subgraph, hits) {
  measures <- data.frame(
    order = vcount(subgraph),
    size = ecount(subgraph),
    density = edge_density(subgraph),
    components = count_components(subgraph),
    Clustering_coefficient = transitivity(subgraph)
  )
  
  degree_values <- degree(subgraph)
  top_deg_index <- order(degree_values, decreasing = TRUE)[1:10]
  top_deg_proteins <- hits$Protein[hits$STRING_id %in% V(subgraph)$name[top_deg_index]]
  
  betweenness_values <- betweenness(subgraph, directed = TRUE, weights = NA)
  top_betweenness_index <- order(betweenness_values, decreasing = TRUE)[1:10]
  top_betweenness_proteins <- hits$Protein[hits$STRING_id %in% V(subgraph)$name[top_betweenness_index]]
  
  eigen_values <- eigen_centrality(subgraph, directed = TRUE, weights = NA)$vector
  top_eigen_index <- order(eigen_values, decreasing = TRUE)[1:10]
  top_eigen_proteins <- hits$Protein[hits$STRING_id %in% V(subgraph)$name[top_eigen_index]]
  
  closeness_values <- closeness(subgraph, mode = "all")
  closeness_values[is.infinite(closeness_values)] <- 0
  top_closeness_index <- order(closeness_values, decreasing = TRUE)[1:10]
  top_closeness_proteins <- hits$Protein[hits$STRING_id %in% V(subgraph)$name[top_closeness_index]]
  
  graph_analysis <- data.frame(
    measures,
    top_deg_proteins = paste(top_deg_proteins, collapse = ";"),
    top_betweenness_proteins = paste(top_betweenness_proteins, collapse = ";"),
    top_eigen_proteins = paste(top_eigen_proteins, collapse = ";"),
    top_closeness_proteins = paste(top_closeness_proteins, collapse = ";")
  )
  
  return(graph_analysis)
}