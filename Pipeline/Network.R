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

igraph_analysis <- function(interactions){
  #Se coge el conjunto de proteínas sobre-expresadas e infra-expresadas devueltas por la función interactions
  hits_upregulated <- as.data.frame(interactions[1])
  hits_downregulated <- as.data.frame(interactions[2])
  
  #Cogemos los 100 identificadores más significativos para formar nuestros subgrafos
  subgraph_up_proteins <- string_db$get_subnetwork(hits_upregulated$STRING_id[1:100])
  subgraph_down_proteins <- string_db$get_subnetwork( hits_downregulated$STRING_id[1:100])
  
  #Representamos dichos subgrafos
  filename = "Interaction graph of upregulated proteins.tiff"
  tiff(filename = filename, units="in", width=9, height=9, res=300)
  plot(subgraph_up_proteins, edge.arrow.size=0.5,vertex.color="gold", vertex.size=5, 
       vertex.frame.color="gray", vertex.label.color="black", vertex.label = hits_upregulated$Protein, 
       vertex.label.cex=.5, vertex.label.dist=2, edge.curved=0.5)
  dev.off()
  filename = "Interaction graph of downregulated proteins.tiff"
  tiff(filename = filename, units="in", width=9, height=9, res=300)
  plot(subgraph_down_proteins, edge.arrow.size=0.5,vertex.color="gold", vertex.size=5, 
       vertex.frame.color="gray", vertex.label.color="black", vertex.label = hits_downregulated$Protein, 
       vertex.label.cex=.5, vertex.label.dist=2, edge.curved=0.5)
  dev.off()
  #########################################################################################
  #Procedemos al análisis por topología
  #########################################################################################
  #Up- regulated
  orden_up <- vcount(subgraph_up_proteins)
  tamaño_up <- ecount(subgraph_up_proteins)
  densidad_up <- edge_density(subgraph_up_proteins) 
  componentes_conexas_up <- count_components(subgraph_up_proteins)
  Clustering_coefficient_up <- transitivity(subgraph_up_proteins)
  
  #Estructuras de parametros grafos
  deg_up <- degree(subgraph_up_proteins)
  top_deg_up <- deg_up[order(deg_up, decreasing = TRUE)[1:10]]
  
  #hist(deg_up, breaks =0:max(deg_up), main = "Histogram of nodes degree", col = "maroon", xlab = "degree")
  #deg.dist <- degree.distribution(subgraph_up_proteins, cumulative = T, mode = "all")
  #plot( x=0:max(deg_up), y=1-deg.dist, pch=19, cex=1.2, col="orange", xlab="Degree", ylab="Cumulative Frequency", main = "Degree distribution")
  
  betweenness_up <- betweenness(subgraph_up_proteins, directed=T, weights=NA)
  top_betweenness_up <- betweenness_up[order(betweenness_up, decreasing = TRUE)[1:10]]
  
  eigen_centrality_up <- eigen_centrality(subgraph_up_proteins, directed=T, weights=NA)
  eigen_centrality_up <- format(eigen_centrality_up$vector, scientific = F, nsmall = 7)
  top_eigen_centrality_up <- order(eigen_centrality_up, decreasing = TRUE)[1:10]
  
  closeness_up <- closeness(subgraph_up_proteins, mode = "all")
  closeness_up <- sub("NaN", 0, closeness_up)
  top_closeness_up <- closeness_up[order(closeness_up, decreasing = TRUE)[1:10]]
  
  #Creamos vectores
  
  nodos_mas_grado_up <- c()
  nodo_mas_eigen_up <- c()
  nodo_mas_betweenness_up <- c()
  nodo_mas_closeness_up <- c()
  
  #iteramos
  for (i in seq(1:10)){
    new_idgrado_up <- as.data.frame(interactions[1])$Protein[as.data.frame(interactions[1])$STRING_id == V(subgraph_up_proteins)$name[degree(subgraph_up_proteins)==top_deg_up[i]]]
    nodos_mas_grado_up <- c(nodos_mas_grado_up, new_idgrado_up)
    new_idbetweenness_up <- as.data.frame(interactions[1])$Protein[as.data.frame(interactions[1])$STRING_id ==V(subgraph_up_proteins)$name[betweenness(subgraph_up_proteins)==top_betweenness_up[i]]]
    nodo_mas_betweenness_up <- c(nodo_mas_betweenness_up, new_idbetweenness_up)
    new_ideigen_up <- as.data.frame(interactions[1])$Protein[as.data.frame(interactions[1])$STRING_id ==V(subgraph_up_proteins)$name[eigen_centrality_up==eigen_centrality_up[top_eigen_centrality_up[i]]]]
    nodo_mas_eigen_up <- c(nodo_mas_eigen_up, new_ideigen_up)
    new_idcloseness_up <- as.data.frame(interactions[1])$Protein[as.data.frame(interactions[1])$STRING_id ==V(subgraph_up_proteins)$name[closeness_up==top_closeness_up[i]]]
    nodo_mas_closeness_up <- c(nodo_mas_closeness_up, new_idcloseness_up)
  }
  nodos_mas_grado_up <- paste(nodos_mas_grado_up, collapse =  ";")
  nodo_mas_betweenness_up <- paste(nodo_mas_betweenness_up, collapse =  ";")
  nodo_mas_eigen_up <- paste(nodo_mas_eigen_up, collapse =  ";")
  nodo_mas_closeness_up <- paste(nodo_mas_closeness_up, collapse =  ";")
  
  graph_analysis_up <- data.frame(orden_up, tamaño_up, densidad_up, componentes_conexas_up, nodos_mas_grado_up, nodo_mas_betweenness_up, nodo_mas_eigen_up, nodo_mas_closeness_up, Clustering_coefficient_up)
  
  
  #Down- regulated
  orden_down <- vcount(subgraph_down_proteins)
  tamaño_down <- ecount(subgraph_down_proteins)
  densidad_down <- edge_density(subgraph_down_proteins) 
  componentes_conexas_down <- count_components(subgraph_down_proteins)
  Clustering_coefficient_down <- transitivity(subgraph_down_proteins)
  
  #Estructuras de parametros grafos
  deg_down <- degree(subgraph_down_proteins)
  top_deg_down <- deg_down[order(deg_down, decreasing = TRUE)[1:10]]
  
  #hist(deg_down, breaks =0:max(deg_down), main = "Histogram of nodes degree", col = "maroon", xlab = "degree")
  #deg.dist <- degree.distribution(subgraph_down_proteins, cumulative = T, mode = "all")
  #plot( x=0:max(deg_down), y=1-deg.dist, pch=19, cex=1.2, col="orange", xlab="Degree", ylab="Cumulative Frequency", main = "Degree distribution")
  
  betweenness_down <- betweenness(subgraph_down_proteins, directed=T, weights=NA)
  top_betweenness_down <- betweenness_down[order(betweenness_down, decreasing = TRUE)[1:10]]
  
  eigen_centrality_down <- eigen_centrality(subgraph_down_proteins, directed=T, weights=NA)
  eigen_centrality_down <- format(eigen_centrality_down$vector, scientific = F, nsmall = 7)
  top_eigen_centrality_down <- order(eigen_centrality_down, decreasing = TRUE)[1:10]
  
  closeness_down <- closeness(subgraph_down_proteins, mode = "all")
  closeness_down <- sub("NaN", 0, closeness_down)
  top_closeness_down <- closeness_down[order(closeness_down, decreasing = TRUE)[1:10]]
  
  #Creamos vectores
  
  nodos_mas_grado_down <- c()
  nodo_mas_eigen_down <- c()
  nodo_mas_betweenness_down <- c()
  nodo_mas_closeness_down <- c()
  
  #iteramos
  for (i in seq(1:10)){
    new_idgrado_down <- as.data.frame(interactions[2])$Protein[as.data.frame(interactions[2])$STRING_id == V(subgraph_down_proteins)$name[degree(subgraph_down_proteins)==top_deg_down[i]]]
    nodos_mas_grado_down <- c(nodos_mas_grado_down, new_idgrado_down)
    new_idbetweenness_down <- as.data.frame(interactions[2])$Protein[as.data.frame(interactions[2])$STRING_id ==V(subgraph_down_proteins)$name[betweenness(subgraph_down_proteins)==top_betweenness_down[i]]]
    nodo_mas_betweenness_down <- c(nodo_mas_betweenness_down, new_idbetweenness_down)
    new_ideigen_down <- as.data.frame(interactions[2])$Protein[as.data.frame(interactions[2])$STRING_id ==V(subgraph_down_proteins)$name[eigen_centrality_down==eigen_centrality_down[top_eigen_centrality_down[i]]]]
    nodo_mas_eigen_down <- c(nodo_mas_eigen_down, new_ideigen_down)
    new_idcloseness_down <- as.data.frame(interactions[2])$Protein[as.data.frame(interactions[2])$STRING_id ==V(subgraph_down_proteins)$name[closeness_down==top_closeness_down[i]]]
    nodo_mas_closeness_down <- c(nodo_mas_closeness_down, new_idcloseness_down)
  }
  nodos_mas_grado_down <- paste(nodos_mas_grado_down, collapse =  ";")
  nodo_mas_betweenness_down <- paste(nodo_mas_betweenness_down, collapse =  ";")
  nodo_mas_eigen_down <- paste(nodo_mas_eigen_down, collapse =  ";")
  nodo_mas_closeness_down <- paste(nodo_mas_closeness_down, collapse =  ";")
  
  graph_analysis_down <- data.frame(orden_down, tamaño_down, densidad_down, componentes_conexas_down, nodos_mas_grado_down, nodo_mas_betweenness_down, nodo_mas_eigen_down, nodo_mas_closeness_down, Clustering_coefficient_down)
  
  graph_analysis <- list(graph_analysis_up, graph_analysis_down)
  return(graph_analysis)
}