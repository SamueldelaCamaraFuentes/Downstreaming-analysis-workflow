Goterms_finder <- function(df, target, numeric_ns, mthreshold, filter_na, organismo, ...){
  
  up <- bind_rows(
    df %>% 
      filter(df$expression == 'Up-regulated') %>% 
      arrange(p.mod) 
  )
  up <- subset(up, up$p.mod <= 0.05)
  
  down <- bind_rows(
    df %>% 
      filter(df$expression == 'Down-regulated') %>% 
      arrange(p.mod) 
  )
  down <- subset(down, down$p.mod <= 0.05)
  
  #Afianzamos la obtención de los identificadores de nuestras proteínas a ensemblegenome
  up_names <- gconvert(up$Protein, organism = "calbicans",  target, numeric_ns, mthreshold, filter_na)
  down_names <- gconvert(down$Protein, organism = "calbicans",  target, numeric_ns, mthreshold, filter_na)
  
  #Multi enrichment analysis
  
  multi_gp <- gost(list("up-regulated" = up_names$name, "down-regulated" = down_names$name), organism = organismo, ...)
  multi_gp[[1]]$adj.P.Val <- p.adjust(multi_gp[[1]]$p_value, method = "fdr")
  # modify the g:Profiler data frame
  gp_mod <- multi_gp$result[, c("query", "source", "term_id",   #Se extraen las columnas de interés
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
  # dedinimos objeto enrichResult
  gp_mod_enrich <- new("enrichResult", result = gp_mod2) 
  
  go_structures <- list(gp_mod_cluster, gp_mod_enrich, multi_gp)
  
  return(go_structures)
}

dotplot_func <- function(Go_terms, ...){
  
  filename = "Dotplot of Go terms.tiff"
  tiff(filename = filename, units="in", width=9, height=9, res=300)
  dotplot(Go_terms[[1]],  ...)
  dev.off()
}

gostplot_func <- function(terms, ...){
  p <- gostplot(terms[[3]], ...) 
  return(p)
  
}

barplot_func <- function(terms, conditions, ...){
  barplot(terms[[2]], ...) + 
    ggplot2::facet_grid(~Cluster) + ggplot2::ylab("Number of proteins") +ggtitle(conditions)
  
  ggsave("barplot.tiff", units="in", width=7, height=7, dpi=300)
}