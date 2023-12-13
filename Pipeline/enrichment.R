Goterms_finder <- function(df, target, numeric_ns, mthreshold, filter_na, organismo, ...){
    #@ df es el dataframe que es la salida de limma.
  up <- bind_rows(
    df %>% 
      filter(df$expression == 'Up-regulated') %>% 
      arrange(p.value) 
  )
  #up <- subset(up, up$p.value <= 0.05)
  
  down <- bind_rows(
    df %>% 
      filter(df$expression == 'Down-regulated') %>% 
      arrange(p.value) 
  )
  #down <- subset(down, down$p.value <= 0.05)
  
  #Afianzamos la obtenci?n de los identificadores de nuestras prote?nas a ensemblegenome
  up_names <- gconvert(up$Protein, organism = organismo,  target, numeric_ns, mthreshold, filter_na)
  down_names <- gconvert(down$Protein, organism = organismo,  target, numeric_ns, mthreshold, filter_na)
  
  #Multi enrichment analysis
  
  multi_gp <- gost(list("up-regulated" = up_names$name, "down-regulated" = down_names$name), organism = organismo, ...)
  multi_gp[[1]]$adj.P.Val <- p.adjust(multi_gp[[1]]$p_value, method = "fdr")
  gp_mod <- multi_gp$result[, c("query", "source", "term_id",   #Se extraen las columnas de inter?s
                                "term_name", "p_value","adj.P.Val", "query_size",
                                "intersection_size" , "term_size",
                                "effective_domain_size", "intersection")]
  gp_mod$ProteinRatio <- paste0(gp_mod$intersection_size, "/", gp_mod$query_size) #Creamos la columna GeneRatio
  
  gp_mod$BgRatio <- paste0(gp_mod$term_size, "/", gp_mod$effective_domain_size)#Creamos la columna BgRatio
  names(gp_mod) <- c("Cluster", "Category", "ID", "Description", "p.value", "adj.P.Val",
                     "query_size", "Count", "term_size", "effective_domain_size",
                     "geneID", "GeneRatio", "BgRatio")
  
  gp_mod$geneID <- gsub(",", "/", gp_mod$geneID)
  gp_mod$Conditions <- case_when( gp_mod$Cluster == "up-regulated" ~ "Up-regulated",
                                  gp_mod$Cluster == "down-regulated" ~ "Down-regulated")
  gp_mod2 <- gp_mod[!duplicated(gp_mod$ID), ]
  row.names(gp_mod2) <- gp_mod2$ID
  
  # definimos objeto compareCluster
  gp_mod_cluster <- new("compareClusterResult", compareClusterResult = gp_mod2)
  # definimos objeto enrichResult
  gp_mod_enrich <- new("enrichResult", result = gp_mod2) 
  
  go_structures <- list(gp_mod_cluster, gp_mod_enrich, multi_gp)
  
  return(go_structures)
}

dotplot_func <- function(Go_terms, ...){
  
  dotplot(Go_terms[[1]],  ...) +  facet_grid(.~Conditions)
}

gostplot_func <- function(terms, ...){
  p <- gostplot(terms[[3]], ...) 
  return(p)
  
}

barplot_func <- function(terms, number, conditions, ...){
  
  
  terms[[2]]@result <- terms[[2]]@result[order(terms[[2]]@result$Count, decreasing = TRUE), ]
  
  #terms_display <- number * 2
  
  
  up_regulated_subset <- terms[[2]]@result[terms[[2]]@result$Conditions == "Up-regulated", ][1:number,]
  down_regulated_subset <- terms[[2]]@result[terms[[2]]@result$Conditions == "Down-regulated", ][1:number,]
  terms[[2]]@result <- rbind(up_regulated_subset, down_regulated_subset)
  
  
  barplot(terms[[2]], showCategory = 100, ...) + 
    ggplot2::facet_grid(~Cluster) + ggplot2::ylab("Number of proteins") +ggtitle(conditions)
  
}