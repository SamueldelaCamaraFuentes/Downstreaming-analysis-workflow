#################################### Pipeline example #################################### 


# Installing CRAN packages:
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

#Installing bioconductor packages: 
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
string_db <- STRINGdb$new(version="11.5", species=237561, score_threshold=200, input_directory="")



packagesdir <- readline("Introduzca el directorio donde se encuentran las funciones:")

setwd(packagesdir)

source("dataprocessing.R")
source("quality_metrics.R")
source("statisticalAnalysis.R")
source("overviewfigures.R")
source("enrichment.R")
source("Network.R")

#Source dir
initialDirectory <- readline("Introduzca el directorio donde estan los resultados de MaxQuant:")

setwd(initialDirectory)

raw <- read.delim("proteinGroups.txt", sep = "\t", stringsAsFactors = FALSE, colClasses = "character") 

newDirectory <- readline('Introduzca el directorio donde generar los resultados:') 

setwd(newDirectory)

dir_create <- readline("Introduzca el nombre del directorio donde guardar los plots:")

dir.create(dir_create) 
finalDirectory = paste0(newDirectory,dir_create)
setwd(finalDirectory)

##################################################################################
#Filtrado inicial
replicas_condicion1 <- as.numeric(readline("Introduzca el número de réplicas biológicas de su condicion control:"))
if (eval(is.na(replicas_condicion1))) {
  replicas_condicion1 <- as.numeric(readline("Error message, Introduzca con número el número de réplicas biológicas de su condición control:"))
}
replicas_condicion2 <- as.numeric(readline("Introduzca el número de réplicas biológicas de su condicion tratamiento:"))
if (eval(is.na(replicas_condicion2))) {
  replicas_condicion2 <- as.numeric(readline("Error message, Introduzca con número el número de réplicas biológicas de su condición tratamiento:"))
}

df<- quick_filtering(raw)

##################################################################################
#Obtención Log2 names

LOG2.names <- obtain_LOG.names(df)

##################################################################################
#Exclusive proteins
first_condition <- readline("Introduzca la expresión regular de su primera condicion:") #WT[0-9]
second_condition <- readline("Introduzca la expresión regular de su segunda condicion:") #WT_H2O2
conditions1 <- c(first_condition, second_condition)
condition1_names <- grep(conditions1[1], LOG2.names, value = TRUE)
condition2_names <- grep(conditions1[2], LOG2.names, value = TRUE)
cond.names <- c(condition1_names,condition2_names)

unique_proteins <- obtain_unique_proteins(df, conditions1, LOG2.names)

##################################################################################
#Venn Diagram

venn <- venn_diagram(df, unique_proteins, fill_color = c("blue", "maroon"))

##################################################################################
#Preprocesamiento

min_number <- as.numeric(readline("Introduzca el minimo numero de replicas biológicas en los que encontrar cada proteina:"))
if (eval(is.na(min_number))) {
  min_number <- as.numeric(readline("Error message, Introduzca el minimo numero de replicas biológicas en los que encontrar cada proteina:"))
} else if (min_number > replicas_condicion1){
  min_number <- as.numeric(readline("Error message, Introduzca un número que no sea mayor al número de réplicas:"))
} else if (min_number > replicas_condicion2){
  min_number <- as.numeric(readline("Error message, Introduzca un número que no sea mayor al número de réplicas:"))
}
min_count <- c(min_number,min_number)

#Filtrado

df.F <- filter_valids(df, unique_proteins, conditions1, min_count, at_least_one <- TRUE, LOG2.names)

#Normalización e imputación

cat("Decida los pasos a proceder:")
cat("1) Normalizar e imputar.")
cat("2) Solo imputar.")
choice <- as.numeric(readLines(con = stdin(), n = 1))
if (eval(is.na(choice))) {
  choice <- as.numeric(readline("Error message, Introduzca un número:"))
}
if (choice == 1){
  df.F <- median_centering(df.F, LOG2.names)
  cat("Decida un método para imputa:")
  cat("1) Imputación distribución normal.")
  cat("2) Imputación por método K-Nearest Neighbours.")
  imp_choice <- as.numeric(readLines(con = stdin(), n = 1))
  if (eval(is.na(imp_choice))) {
    imp_choice <- as.numeric(readline("Error message, Introduzca un número:"))
  }
  if (imp_choice == 1){
    df.FNI <- impute_data(df.F, LOG2.names)
  } else if (imp_choice == 2){
    df.FNI <- impute_KNN_data(df.F, LOG2.names, k = 5)
  }
} else if (choice == 2){
  cat("Decida un método para imputar:")
  cat("1) Imputación distribución normal.")
  cat("2) Imputación por método K-Nearest Neighbours.")
  imp_choice2 <- as.numeric(readLines(con = stdin(), n = 1))
  if (imp_choice2 == 1){
    df.FNI <- impute_data(df.F, LOG2.names)
  } else if (imp_choice2 == 2){
    df.FNI <- impute_KNN_data(df.F, LOG2.names, k = 5)
  }
  
}

#Control

plotCV2(df.FNI[,LOG2.names],  trend = TRUE, main = "Dispersion check", cex = 0.2, pch = 16, xlab="Average log-intensity", ylab=expression("Relative standard deviation"))

boxplot <- boxplot_function(df.FNI, cond.names, cex.axis = 0.5)

imputation_state(df.F, df.FNI, cond.names)

histogram(df.FNI, cond.names, color = "maroon")

setwd(finalDirectory)


##################################################################################
#Expresión diferencial

limma <- limma.analysis(df.FNI) #LOG2.WT1, LOG2.WT_H2O2_1, LOG2.Prn1_1 LOG2.Prn1_H2O2_1

##################################################################################
#Representaciones
conditions <- readline("Introduzca el titulo de sus condiciones a enfrentar:") #"prn1\u2206 vs prn1\u2206 H2O2 10mM"

volcano <- volcano_plot(limma, conditions) 

my_pca <- pca(df.FNI,cond.names)

top_proteins <- my_heatmap(limma, df.FNI, cond.names, conditions, whole = TRUE, Rowv = NULL, Colv = NA, col =greenred(75), cexCol = 0.6 )

##################################################################################
#Análisis funcional
id_target <- readline("Introducir el tipo de identificador proteico al que mapear:") #ENSG
organismo <- readline("Introducir organismo:") #calbicans
threshold <- as.numeric(readline("Introduzca un valor umbral para la significancia:"))
if (eval(is.na(threshold))) {
  threshold <- as.numeric(readline("Error message, Introduzca un valor umbral numérico para la significancia:"))
}
Go_terms <- Goterms_finder(limma, target = id_target, numeric_ns = "", mthreshold = Inf, filter_na = TRUE, organism = organismo, user_threshold = threshold, multi_query = FALSE, evcodes = TRUE, sources = c("GO", "KEGG", "WP", "REAC", "CORUM"))

font_size <- as.numeric(readline("Introduzca el tamaño fuente:"))
if (eval(is.na(font_size))) {
  font_size <- as.numeric(readline("Error message, Introduzca un tamaño fuente numérico:"))
}
terms_number <- as.numeric(readline("Introduzca el número de términos a mostrar:"))
if (eval(is.na(terms_number))) {
  terms_number <- as.numeric(readline("Error message, Introduzca un número de términos a mostrar numérico:"))
}
dotplot_func(Go_terms, x = "Conditions", title = conditions, font.size = font_size, showCategory = terms_number, color = "adj.P.Val")

p <- gostplot_func(Go_terms, interactive = FALSE)

#Descarga del gráfico de manhattan
publish_gostplot(p, filename = "manhattan.tiff")

#Descarga de la tabla con los términos
publish_gosttable(Go_terms[[3]]$result, filename = "goterms.tiff")

barplot_func(Go_terms, conditions, showCategory = 100, font.size = 10)

##################################################################################
#Análisis de interacción de proteínas

interactions_analysis <- interactions(limma)

graph_analysis <- igraph_analysis(interactions_analysis)
