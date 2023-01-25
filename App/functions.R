
# Installing CRAN packages:
if(!require(BiocManager)){install.packages("BiocManager")}
library(BiocManager)
if(!require(methods)){install.packages("methods")}
library(methods)
if(!require(colorspace)){install.packages("colorspace")}
library(colorspace)
if(!require(grid)){install.packages("grid")}
library(grid)
if(!require(utils)){install.packages("utils")}
library(utils)
if(!require(dplyr)){install.packages("dplyr")}
library(dplyr)
if(!require(ggplot2)){install.packages("ggplot2")}
library(ggplot2)
if(!require(ggvenn)){install.packages("ggvenn")}
library(ggvenn)
if(!require(VIM)){install.packages("VIM")}
library(VIM)
if(!require(gplots)){install.packages("gplots")}
library(gplots)
if(!require(gprofiler2)){install.packages("gprofiler2")}
library(gprofiler2) 
if(!require(writexl)){install.packages("writexl")}
library(writexl)
if(!require(igraph)){install.packages("igraph")}
library(igraph)
if(!require(plotly)){install.packages("plotly")}
library(plotly)

#Installing bioconductor packages: 
if(!require(rsconnect)){BiocManager::install("rsconnect",update=F,ask=F)}
library(rsconnect)
if(!require(BiocGenerics)){BiocManager::install("BiocGenerics",update=F,ask=F)}
library(BiocGenerics)
if(!require(Biobase)){BiocManager::install("Biobase",update=F,ask=F)}
library(Biobase)
if(!require(S4Vectors)){BiocManager::install("S4Vectors",update=F,ask=F)}
library(S4Vectors)
if(!require(IRanges)){BiocManager::install("IRanges",update=F,ask=F)}
if(!require(AnnotationDbi)){BiocManager::install("AnnotationDbi",update=F,ask=F)}
library(AnnotationDbi)
if(!require(limma)){BiocManager::install("limma",update=F,ask=F)}
library(limma)
if(!require(qvalue)){BiocManager::install("qvalue",update=F,ask=F)}
library(qvalue)
if(!require(clusterProfiler)){BiocManager::install("clusterProfiler",update=F,ask=F)}
library(clusterProfiler)
if(!require(enrichplot)){BiocManager::install("enrichplot",update=F,ask=F)}
library(enrichplot)
if(!require(DOSE)){BiocManager::install("DOSE",update=F,ask=F)}
library(DOSE)
if(!require(STRINGdb)){BiocManager::install("STRINGdb",update=F,ask=F)}
library(STRINGdb)

string_db <- STRINGdb$new(version="11.5", species=237561, score_threshold=400, input_directory="", protocol="http")


#################################### helper functions #################################### 


##################################################################################
#Quick filtering
quick_filtering <- function(raw, input){
  df <- raw %>%
    filter(Potential.contaminant != "+") %>% #Nos quedamos con todos aquellos que no tengan un +
    filter(Reverse != "+") %>% #Nos quedamos con todos aquellos que no tengan un +
    filter(Only.identified.by.site != "+") #Nos quedamos con todos aquellos que no tengan un +
  
  #Obtenemos los identificadores
  regex <- regexpr("(?<=orf19.).*(?=.CGDID)", df$Fasta.headers, perl = TRUE)
  df$Protein <- regmatches(df$Fasta.headers, regex)
  
  df$Protein <- sub(" ", ";", df$Protein)
  regex2 <- regexpr("(?<=;).*(?>[A-Z0-9])", df$Protein, perl = TRUE)
  df$Protein <- regmatches(df$Protein, regex2)
  
  #Obtenemos la descripción proteica
  regex3 <- regexpr("(?<=;).*(?>;|[a-z])", df$Fasta.headers, perl = TRUE)
  df$Protein_description <- regmatches(df$Fasta.headers, regex3)
  
  #Hacemos la log-transformación
  if (input == "lfq"){
    intensity_names <- grep("LFQ.intensity", colnames(df), value = TRUE)
    df[intensity_names] <- sapply(df[intensity_names], as.numeric) #convertimos los valores en numericos
    LOG2.names <- sub("^LFQ.intensity", "LOG2", intensity_names)
    df[LOG2.names] <- log2(df[intensity_names])
  } else if (input == "int"){
    intensity_names <- grep("^Intensity", colnames(df), value = TRUE)
    df[intensity_names] <- sapply(df[intensity_names], as.numeric) #convertimos los valores en integers
    LOG_names <- sub("^Intensity", "LOG2", intensity_names)
    df[LOG_names] <- log2(df[intensity_names])
    
  }
  
  return(df)
}

##################################################################################
#Obtaining LOG2 names

obtain_LOG.names <- function(df){
  
  LOG2.names <- grep("LOG2", colnames(df))
  return(colnames(df[LOG2.names]))
}

##################################################################################
#Subseting unique proteins for each condition

obtain_unique_proteins <- function(df, conditions, LOG2.names, replicas_condicion1, replicas_condicion2){
  cond.names <- lapply(conditions, # Sobre la lista conditions, aplicamos una funcion quedarnos con las columnas de las condiciones con datos en log
                       function(x) grep(x, LOG2.names, value = TRUE, perl = TRUE))
  condi_names <- c(cond.names[[1]], cond.names[[2]])
  
  df2 <- df[condi_names]   # Extrae las columnas de interés
  df2 <- as.matrix(df2)
  finite_sums <- rowSums(is.finite(df2[, condi_names[1:replicas_condicion1]])) # Cuenta el numero de valores validos para cada condicion 
  infinite_sums <- rowSums(!is.finite(df2[,condi_names[(replicas_condicion1+1):(replicas_condicion1+replicas_condicion2)]]))
  df$cond1_exclusive <- case_when(finite_sums == replicas_condicion1 & infinite_sums == replicas_condicion2 ~ TRUE,
                                  finite_sums == (replicas_condicion1-1) & infinite_sums == replicas_condicion2 ~ TRUE,
                                  finite_sums== (replicas_condicion1-2) & infinite_sums == replicas_condicion2 ~ TRUE,
                                  TRUE ~ FALSE)
  
  df$cond2_exclusive <- case_when(finite_sums == 0 & infinite_sums == 0 ~ TRUE,
                                  finite_sums == 0 & infinite_sums == replicas_condicion2-(replicas_condicion2-1) ~ TRUE,
                                  finite_sums == 0 & infinite_sums == replicas_condicion2-(replicas_condicion2-2) ~ TRUE,
                                  TRUE ~ FALSE)
  
  
  cond1_unicas <- filter(df, df$cond1_exclusive)
  cond2_unicas <- filter(df, df$cond2_exclusive)
  
  return(list(cond1_unicas, cond2_unicas))
}

##################################################################################
#Venn Diagram

venn_diagram <- function(df, unique_proteins, label1, label2, ...){
  #Proteinas en muestras WT
  cond1_set <- df%>%
    filter(df$Protein != unique_proteins[[2]]$Protein[1])
  for (i in unique_proteins[[2]]$Protein){
    cond1_set <- cond1_set%>%
      filter(cond1_set$Protein != i)
  }
  
  #Proteinas en muestras WT_H2O2
  cond2_set <- df%>%
    filter(df$Protein != unique_proteins[[1]]$Protein[1])
  for (i in unique_proteins[[1]]$Protein){
    cond2_set <- cond2_set %>%
      filter(cond2_set$Protein != i)
  }
  
  x <- list( cond1_set$Protein,  cond2_set$Protein)
  names(x) <- c(label1, label2)
  ggvenn(x, ...)
}


##################################################################################
#Functional analysis of exclusive proteins

exclusive_proteins_GO <- function(unique_proteins, target, numeric_ns, mthreshold, filter_na, ...){
  cond1_names <- gconvert(unique_proteins[[1]]$Protein, organism = "calbicans", target, numeric_ns, mthreshold, filter_na)
  cond2_names <- gconvert(unique_proteins[[2]]$Protein, organism = "calbicans", target, numeric_ns, mthreshold, filter_na)
  
  #Multi enrichment analysis
  
  multi_gp <- gost(list("condition1" = cond1_names$name, "condition2" = cond2_names$name), ...)
  multi_gp[[1]]$adj.P.Val <- p.adjust(multi_gp[[1]]$p_value, method = "fdr")
  multi_gp[[1]] <- arrange(multi_gp[[1]], multi_gp[[1]]$p_value)[1:50,]
  
  # modify the g:Profiler data frame
  gp_mod <- multi_gp$result[, c("query", "source", "term_id",   #Se extraen las columnas de interés
                                "term_name", "p_value","adj.P.Val", "query_size",
                                "intersection_size" , "term_size",
                                "effective_domain_size", "intersection")]
  gp_mod$ProteinRatio <- paste0(gp_mod$intersection_size, "/", gp_mod$query_size) #Creamos la columna GeneRatio
  
  gp_mod$BgRatio <- paste0(gp_mod$term_size, "/", gp_mod$effective_domain_size)#Creamos la columna BgRatio
  names(gp_mod) <- c("Cluster", "Category", "ID", "Description", "p.adjust", "adj.P.Val",
                     "query_size", "Count", "term_size", "effective_domain_size",
                     "geneID", "GeneRatio", "BgRatio")
  
  gp_mod$geneID <- gsub(",", "/", gp_mod$geneID)
  gp_mod$Conditions <- case_when( gp_mod$Cluster == "up-regulated" ~ "Up-regulated",
                                  gp_mod$Cluster == "down-regulated" ~ "Down-regulated")
  row.names(gp_mod) <- make.names(gp_mod$ID, unique = TRUE) #Utilizamos make.names para el caso de que se produzca una redundancai de identificadores
  
  # definimos objeto compareCluster
  gp_mod_cluster <- new("compareClusterResult", compareClusterResult = gp_mod)
  # dedinimos objeto enrichResult
  gp_mod_enrich <- new("enrichResult", result = gp_mod) 
  
  go_structures <- list(gp_mod_cluster, gp_mod_enrich, multi_gp)
  write_xlsx(gp_mod, paste0(newDirectory,"/plots/Analisis_funcional_proteinas_exclusivas.xlsx"))
  return(go_structures)
}


##################################################################################
#Filtering

filter_valids <- function(df, unique_proteins, conditions1, min_count, at_least_one = FALSE, LOG2.names) {
  
  
  cond.names <- lapply(conditions1, # Sobre la lista conditions, aplicamos una funcion quedarnos con las columnas de las condiciones con datos en log
                       function(x) grep(x, LOG2.names, value = TRUE, perl = TRUE))
  print(cond.names)
  
  total_unique_proteins <- rbind(unique_proteins[[1]], unique_proteins[[2]])
  common_df <- df%>%
    filter(df$Protein != total_unique_proteins$Protein[1])
  for (i in total_unique_proteins$Protein){
    common_df <- common_df%>%
      filter(common_df$Protein != i)
  }
  
  cond.filter <- sapply(1:length(cond.names), function(i) { #En nuestro caso se itera del 1 al 2
    df2 <- common_df[cond.names[[i]]]   # Extrae las columnas de interés
    df2 <- as.matrix(df2)   # Lo convierte en matriz
    sums <- rowSums(is.finite(df2)) # Cuenta el numero de valores validos para cada condicion 
    sums >= min_count[i]   # Calculates whether min_count requirement is met si el número de datos es igual o superior al elemento en min_count se establece true
  })
  if (at_least_one) {
    common_df$KEEP <- apply(cond.filter, 1, any)
  } else {
    common_df$KEEP <- apply(cond.filter, 1, all)
  }
  common_df <- filter(common_df, KEEP)
  
  common_df[LOG2.names] <- lapply(LOG2.names,
                           function(x) {
                             temp <- common_df[[x]]
                             temp[!is.finite(temp)] = NA
                             return(temp)
                             
                           })
  
  
  
  return(common_df)  
}

##################################################################################
#Normalization 

median_centering <- function(df, LOG2.names) {
  
  df[, LOG2.names] <- lapply(LOG2.names, #A las columnas de interés, se procede a calcular el LOG2-mediana de cada 
                             function(x) {
                               LOG2 <- df[[x]] #Se asigna el valor calculado
                               LOG2[!is.finite(LOG2)] <- NA   # Excluimos los missing values del cálculo de la mediana 
                               gMedian <- median(LOG2, na.rm = TRUE)
                               LOG2 - gMedian #Se calcula el valor.
                             }
  )
  
  return(df)
}



##################################################################################
#Imputation

impute_data <- function(df, LOG2.names, width = 0.3, downshift = 1.8) {
  
  impute.names <- sub("^LOG2", "impute", LOG2.names)
  
  # Create new columns indicating whether the values are imputed
  df[impute.names] <- lapply(LOG2.names, function(x) !is.finite(df[, x]))
  
  # Imputation
  set.seed(1)
  df[LOG2.names] <- lapply(LOG2.names,
                           function(x) {
                             temp <- df[[x]]
                             temp[!is.finite(temp)] = NA
                             
                             temp.sd <- width * sd(temp[df$KEEP], na.rm = TRUE)   # shrink sd width
                             temp.mean <- mean(temp[df$KEEP], na.rm = TRUE) - 
                               downshift * sd(temp[df$KEEP], na.rm = TRUE)   # shift mean of imputed values
                             
                             n.missing <- sum(is.na(temp))
                             temp[is.na(temp)] <- rnorm(n.missing, mean = temp.mean, sd = temp.sd)                          
                             return(temp)
                           })
  
  
  return(df)
}

impute_KNN_data <- function(df, LOG2.names, ...){
  
  impute.names <- sub("^LOG2", "impute", LOG2.names)
  df[impute.names] <- lapply(LOG2.names, function(x) !is.finite(df[, x])) #Creamos columnas donde se ve si imputar o no
  
  #Imputación
  df[LOG2.names] <- lapply(LOG2.names,
                           function(x) {
                             temp <- df[[x]]
                             temp[!is.finite(temp)] = NA
                             return(temp)
                             
                           })
  
  imp_knn <- kNN(df, variable = LOG2.names, ... )

  return(imp_knn)
}


##################################################################################
#Processing check

plotCV2 <- function(y, trend = TRUE, main= "Imputation check", ...){

  A <- rowMeans(y, na.rm = TRUE) #Hago las medias de las filas
  CV <- (matrixStats::rowSds(data.matrix(y), na.rm = TRUE)/A)^2 #Calculo de la desviación estandar relativa invocando el paquete matrixStats y usando su función rowSds
  res <- data.frame(mean = A, CV = CV)
  plot(A, CV ,  ylim = c(min(CV)-0.001, max(CV) +0.001),  ...)
  if(trend){ 
    fit <- limma::loessFit(CV, A) #Invocamos la funcion loess del paquete limma, para hacer una regresión local
    o <- order(A)
    lines(A[o], fit$fitted[o], lwd =2, col = "red")
  }
  
  return(res)
}

boxplot_function <- function(df, cond.names,  ...){
  par(mfrow=c(1,1), font.lab=2, cex.lab=1, font.axis=2, cex.axis=1, cex.main=1, las = 1)
  boxplot(df[, cond.names], names = cond.names,  ylim = c(min(df[,cond.names])-1, max(df[,cond.names]) +1), main="Boxplot normalized Intensities", ylab = "Intensidades", las=2,  ...)
}

preimputation_state <- function(df, cond.names){

  par(mfrow=c(2,2))
  aggr(df[, cond.names], delimiter = "NA", labels = names(df[cond.names]), cex.axis = .5)
}

postimputation_state <- function(df, cond.names){

  par(mfrow=c(2,2))
  aggr(df[, cond.names], delimiter = "NA", labels = names(df[cond.names]), cex.axis = .5)
}

histogram <- function(df, colname, color, title){
  hist(as.numeric(df[,colname]), main = title, xlab = 'Intensity-value', col = color)
  
  
}

##################################################################################
#Differential expression analysis

limma.analysis <- function(df, LOG2.names, paired = FALSE, replicas_condicion1, replicas_condicion2, condition1, condition2, orden, logfcup, logfcdown, sig, adjval){
  condition1_names <- grep(condition1, LOG2.names, value = TRUE)
  condition2_names <- grep(condition2, LOG2.names, value = TRUE)
  
  #Control columns
  for (i in 1:replicas_condicion1){
    nam <- paste("control", i, sep = "")
    assign(nam, condition1_names[i])
  }  #LOG2.WT1,#LOG2.WT2,#LOG2.WT3,#LOG2.WT4
  
  #Treatment columns
  for (i in 1:replicas_condicion2){
    nam <- paste("prob_column", i, sep = "")
    assign(nam, condition2_names[i])
  }  #LOG2.WT_H2O2_1, #LOG2.WT_H2O2_2, #LOG2.WT_H2O2_3, #LOG2.WT_H2O2_4
  ct <- c()
  for (i in ls()[grep("control", ls())]){
    new_value_control <- get(i)
    ct <- c(ct, new_value_control)
  }
  print(ct)
  
  tr <- c()
  for (i in ls()[grep("prob_column", ls())]){
    new_value_prob <- get(i)
    tr <- c(tr, new_value_prob)
  }
  print(tr)
  
  if (paired == FALSE){
  control <- rep(1, replicas_condicion1)
  treatment <- rep(2, replicas_condicion2)
  design <- model.matrix(~factor(c(control, treatment)))
  } else if (paired == TRUE){
    pairinfo = factor(rep(1:replicas_condicion1,2))
    control <- rep(1, replicas_condicion1)
    treatment <- rep(2, replicas_condicion2)
    design <- model.matrix(~pairinfo+factor(c(control, treatment)))
  }
  
  dat <- df[, c(ct, tr)]
  n <- dim(dat)[1]
  fit <- lmFit(dat, design)
  fit.eb <- eBayes(fit)
  print(colnames(fit.eb))
  logFC <- fit.eb$coefficients[, 2] #Cálculo del log fold-change
  p.mod <- fit.eb$p.value[, 2]    # p-valor moderado correspondiente al estadístico t moderado.
  expression <- case_when(logFC >= logfcup & -log10(p.mod) >= sig ~ "Up-regulated",
                          logFC <= logfcdown & -log10(p.mod) >= sig ~ "Down-regulated",
                          TRUE ~ "Unchanged")#labels para expresion 
  adj.P.Val <- p.adjust(p.mod, method = adjval)
  results.eb <- data.frame(logFC, p.mod, adj.P.Val, expression)
  if (orden == 1){
    results.eb <- results.eb[order(results.eb$p.mod), ] #Ordenamos por p.mod
  } else if (orden ==2){
    results.eb <- results.eb[order(results.eb$adj.P.Val), ]
  }
  results_rownames <-rownames(results.eb) #Obtenemos el nombre de las filas
  Protein <- c()                          #Creamos un vector vacío
  Protein_description <- c()
  for (i in results_rownames){
    new_value <- df[i, "Protein"]     #iteramos los rownames en la columna deseada
    Protein <- c(Protein, new_value)      #Creamos la columna con identificadores
    new_desc <- df[i, "Protein_description"]
    Protein_description <- c(Protein_description, new_desc) #Creamos la columna con descriptores. 
  }
  results.eb$Protein <- Protein 
  results.eb$Protein_description <- Protein_description
  row.names(results.eb) <- results.eb$Protein
  return(results.eb)
}

##################################################################################
#Plots

volcano_plot <- function(limma, title, label){
  
  #top <- readline("Introduzca cuantas proteinas quiere etiquetar:")
  top <- label
  logFC <- limma$logFC
  p.ord <- limma$p.ord
  p.mod <- limma$p.mod
  
  top_proteins <- bind_rows(
    limma %>% 
      filter(limma$expression == 'Up-regulated') %>% #Se cogen las up-regulated y acto seguido se ordenan por p.mod
      arrange(p.mod) %>% 
      head(top),
    limma %>% 
      filter(limma$expression == 'Down-regulated') %>% 
      arrange(p.mod) %>% 
      head(top)
  )
  
  rx <- c(-1, 1)*max(abs(logFC))*1.1 #Para definir nuestro eje X cogemos el máximo de la distribución
  ry <- c(0, ceiling(max(-log10(p.mod), -log10(p.mod)))) #Para definir nuestro eje Y cogemos el el mínimo y el máximo de la distribución
  
  ggplot(limma, aes(logFC, -log10(p.mod))) +
    geom_point(aes(color = expression), size = 2) +
    scale_color_manual(values = c("green3", "grey78", "firebrick3")) +
    guides(colour = guide_legend(override.aes = list(size=1.5))) +
    ggtitle(title) +
    labs(y = "-log10 p-value", x = "log2 Fold change") +
    geom_label(data = top_proteins,
               mapping = aes(logFC, -log10(p.mod), label = Protein),
               size = 1.5)
  
}

pca <- function(x, group){
  kmc<-length(unique(group))
  tdata<-t(x[group])
  km <-kmeans(tdata,(kmc/4),2000) 
  pc<-prcomp(tdata)
  plot(pc$x[,1], pc$x[,2], pch=16, col = c("salmon", "blue"),main="K-means Clustering of PCA", xlab="PC1", ylab="PC2")
  text(pc$x[,1], pc$x[,2],labels=rownames(pc$x),pos=3,offset=0.4,cex=0.7)
  #pc<-cbind(pc$x[,1], pc$x[,2])
  #ordispider(pc, factor(km$cluster), label = TRUE)
  #ordihull(pc, factor(km$cluster), lty = "dotted")
}



my_heatmap <- function(data, cond.names, title){
  
  heatmap <- heatmap(as.matrix(data[cond.names]), labRow = data$Protein, main = title, Rowv = NULL, Colv = NA, col =greenred(75), cexCol = 0.6)
  legend(x="topleft", legend=c("2", "0", "-2"),fill=c("red", "black", "green"), title = "Log2FC")
  return(heatmap)
  
}

my_heatmap_differential <- function(limma, data, cond.names, title){
  top_proteins <- bind_rows(   #Nos quedamos con las proteínas sobre e infra expresadas en las comparativas WT
    limma %>% 
      filter(limma$expression == 'Up-regulated'),  
    limma %>% 
      filter(limma$expression == 'Down-regulated') 
  )
  
  diferential <- data.frame()
  for (i in top_proteins$Protein){   #Nos quedamos con la entrada del dataframe del dataset original de cada proteina de interés
    new <- as.data.frame(data[data$Protein==i,])
    diferential <- rbind(diferential, new)
  }
  diferential[cond.names] <- sapply(diferential[cond.names], as.numeric)
  heatmap <- heatmap(as.matrix(diferential[cond.names]),  labRow = diferential$Protein, main =title, Rowv = NULL, Colv = NA, col =greenred(75), cexCol = 0.6)
  legend(x="left", legend=c("2", "0", "-2"),fill=c("red", "black", "green"), title = "Log2FC")
}

##################################################################################
#Functional analysis



Goterms_finder <- function(df, target, numeric_ns, mthreshold, filter_na, organismo, ...){
  #@ df es el dataframe que es la salida de limma.
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
  up_names <- gconvert(up$Protein, organism = organismo,  target, numeric_ns, mthreshold, filter_na)
  down_names <- gconvert(down$Protein, organism = organismo,  target, numeric_ns, mthreshold, filter_na)
  
  #Multi enrichment analysis
  
  multi_gp <- gost(list("up-regulated" = up_names$name, "down-regulated" = down_names$name), organism = organismo, ...)
  multi_gp[[1]]$adj.P.Val <- p.adjust(multi_gp[[1]]$p_value, method = "fdr")
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
  # definimos objeto enrichResult
  gp_mod_enrich <- new("enrichResult", result = gp_mod2) 
  
  go_structures <- list(gp_mod_cluster, gp_mod_enrich, multi_gp)
  
  return(go_structures)
}

dotplot_func <- function(terms, ...){

  dotplot(terms[[1]],  ...)

}

gostplot_func <- function(terms, ...){
  gostplot(terms[[3]], ...) 
  
}

barplot_func <- function(terms, conditions, ...){

  barplot(terms[[2]], ...) + 
    ggplot2::facet_grid(~Cluster) + ggplot2::ylab("Number of proteins") +ggtitle(conditions)

}

##################################################################################
#Interaction network analysis

interactions_up <- function(df){
  #Subset de up-regulated y down-regulated
  
  up_regulated <- subset(df, df$expression == "Up-regulated")
  down_regulated <- subset(df, df$expression == "Down-regulated")
  
  up_mapped <- string_db$map(up_regulated, "Protein", removeUnmappedRows = TRUE)
  
  down_mapped <- string_db$map(down_regulated, "Protein", removeUnmappedRows = TRUE)
  
  par(mfrow=c(1,1))
  hits_up <- up_mapped$STRING_id[1:100] 
  string_db$plot_network(hits_up)

  interactions <- list(up_mapped, down_mapped)
  return(interactions)
  
}

interactions_down <- function(df){
  down_regulated <- subset(df, df$expression == "Down-regulated")
  
  down_mapped <- string_db$map(down_regulated, "Protein", removeUnmappedRows = TRUE)
  
  par(mfrow=c(1,1))
  hits_down <- down_mapped$STRING_id[1:100] 
  string_db$plot_network(hits_down)


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



