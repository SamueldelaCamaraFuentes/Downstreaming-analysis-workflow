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
if(!require(VennDiagram)){install.packages("VennDiagram")}
library(VennDiagram)
if(!require(VIM)){install.packages("VIM")}
library(VIM)
if(!require(wrProteo)){install.packages("wrProteo")}
library(wrProteo)
if(!require(wrMisc)){install.packages("wrMisc")}
library(wrMisc)
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
if(!require(matrixStats)){install.packages("matrixStats")}
library(matrixStats)
if(!require(ggExtra)){install.packages("ggExtra")}
library(ggExtra)
if(!require(corrplot)){install.packages("corrplot")}
library(corrplot)

#Installing bioconductor packages: 
if(!require(rsconnect)){BiocManager::install("rsconnect",update=F,ask=F)}
library(rsconnect)
if(!require(BiocGenerics)){BiocManager::install("BiocGenerics",update=F,ask=F)}
library(BiocGenerics)
if(!require(Biobase)){BiocManager::install("Biobase",update=F,ask=F)}
library(Biobase)
if(!require(S4Vectors)){BiocManager::install("S4Vectors",update=F,ask=F)}
library(S4Vectors)
if(!require(vsn)){BiocManager::install("vsn",update=F,ask=F)}
library(vsn)
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



#################################### helper functions #################################### 


##################################################################################
#Quick filtering
quick_filtering <- function(raw, input, platform, organism){
  if (platform == 1){
    df <- raw %>%
      filter(Potential.contaminant != "+") %>% #Nos quedamos con todos aquellos que no tengan un +
      filter(Reverse != "+") %>% #Nos quedamos con todos aquellos que no tengan un +
      filter(Only.identified.by.site != "+") #Nos quedamos con todos aquellos que no tengan un +
    
    #Obtenemos los identificadores
    regex <- regexpr(".*(?=.CGDID)", df$Fasta.headers, perl = TRUE)
    df$Protein <- regmatches(df$Fasta.headers, regex)
    
    df$Protein <- sub(" ", ";", df$Protein)
    regex2 <- regexpr("(?<=;).*(?>[A-Z0-9])", df$Protein, perl = TRUE)
    df$Protein <- regmatches(df$Protein, regex2)
    
    #Obtenemos la descripci?n proteica
    regex3 <- regexpr("(?<=;).*(?>;|[a-z])", df$Fasta.headers, perl = TRUE)
    df$Protein_description <- regmatches(df$Fasta.headers, regex3)
    
    #Hacemos la log-transformaci?n
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
      
    } else if (platform == 2){
      df <- raw 
      
      if (organism == 1){
        regex <- regexpr(".*(?=.CGDID)", df$Protein.ID, perl = TRUE)
        df$Protein <- regmatches(df$Protein.ID, regex)
        
        df$Protein <- sub(" ", ";", df$Protein)
        regex2 <- regexpr("(?<=;).*(?>[A-Z0-9])", df$Protein, perl = TRUE)
        df$Protein <- regmatches(df$Protein, regex2)
        
        #Obtenemos la descripci?n proteica
        regex3 <- regexpr("(?<=;).*(?>;|[a-z])", df$Protein.ID, perl = TRUE)
        df$Protein_description <- regmatches(df$Protein.ID, regex3) 
        
        colnames(df)[1] <- "Protein"
        
      } else if (organism == 2){
        
        colnames(df)[1] <- "Protein.ID"
        colnames(df)[2] <- "Protein"
        colnames(df)[3] <- "Protein_description"
      }
  
      if (input == "lfq"){
        intensity_names <- grep("MaxLFQ.Intensity", colnames(df), value = TRUE)
        df[intensity_names] <- sapply(df[intensity_names], as.numeric) #convertimos los valores en numericos
        LOG2.names <- sub("MaxLFQ.Intensity", "LOG2", intensity_names)
        #df[LOG2.names] <- log2(df[intensity_names])
        df[LOG2.names] <- lapply(df[intensity_names], log2)
      } else if (input == "int"){
        intensity_names <- grep("Intensity", colnames(df), value = TRUE)
        df[intensity_names] <- sapply(df[intensity_names], as.numeric) #convertimos los valores en integers
        LOG2.names <- sub("Intensity", "LOG2", intensity_names)
        df[LOG2.names] <- log2(df[intensity_names])
      }
      return(df)
    } else if (platform == 3){
      df <- raw 
      colnames(df)[3] <- "Protein"
      colnames(df)[1] <- "Protein_description"
      
      intensity_columns <- grepl("\\.d", colnames(df)) #Nos quedamos con aquellas columnas que corresponden con muestras
      intensity_columns <- colnames(df[intensity_columns])
      
      #Le damos un lavado 
      new_names <- c()
      for (i in intensity_columns){
        new <- basename(i)
        new_names <- c(new_names, new)
      }
      
      colnames(df)[length(colnames(df)) - (length(intensity_columns)-1):length(colnames(intensity_columns))] <- new_names
      intensity_names <- grep("\\.d", colnames(df), value = TRUE)
      df[intensity_names] <- sapply(df[intensity_names], as.numeric) #convertimos los valores en integers
      LOG2.names <- sub("\\.d", ".LOG2", intensity_names)
      df[LOG2.names] <- log2(df[intensity_names])
      
      
      return(df)
    }
  

}

##################################################################################
#Obtaining LOG2 names

obtain_LOG.names <- function(df){
  
  LOG2.names <- grep("LOG2", colnames(df), value = TRUE)
  return(LOG2.names)
}

##################################################################################
#Subseting unique proteins for each condition

obtain_unique_proteins <- function(df, conditions, LOG2.names, replicas_condicion1, replicas_condicion2){
  cond.names <- lapply(conditions, # Sobre la lista conditions, aplicamos una funcion quedarnos con las columnas de las condiciones con datos en log
                       function(x) grep(x, LOG2.names, value = TRUE, perl = TRUE))
  condi_names <- c(cond.names[[1]], cond.names[[2]])
  
  df2 <- df[condi_names]   # Extrae las columnas de inter?s
  df2 <- as.matrix(df2)
  finite_sums <- rowSums(is.finite(df2[, condi_names[1:replicas_condicion1]])) # Cuenta el numero de valores validos para cada condicion 
  infinite_sums <- rowSums(!is.finite(df2[,condi_names[(replicas_condicion1+1):(replicas_condicion1+replicas_condicion2)]]))
  df$cond1_exclusive <- case_when(finite_sums == replicas_condicion1 & infinite_sums == replicas_condicion2 ~ TRUE,
                                  finite_sums == (replicas_condicion1-1) & infinite_sums == replicas_condicion2 ~ TRUE,
                                  finite_sums== (replicas_condicion1-2) & infinite_sums == replicas_condicion2 ~ TRUE,
                                  finite_sums== (replicas_condicion1-3) & infinite_sums == replicas_condicion2 ~ TRUE,
                                  TRUE ~ FALSE)
  
  df$cond2_exclusive <- case_when(finite_sums == 0 & infinite_sums == 0 ~ TRUE,
                                  finite_sums == 0 & infinite_sums == replicas_condicion2-(replicas_condicion2-1) ~ TRUE,
                                  finite_sums == 0 & infinite_sums == replicas_condicion2-(replicas_condicion2-2) ~ TRUE,
                                  finite_sums == 0 & infinite_sums == replicas_condicion2-(replicas_condicion2-3) ~ TRUE,
                                  TRUE ~ FALSE)
  
  
  cond1_unicas <- filter(df, df$cond1_exclusive)
  cond2_unicas <- filter(df, df$cond2_exclusive)
  
  return(list(cond1_unicas, cond2_unicas))
}

##################################################################################
#Proteins identified
identify_proteins <- function(raw, condi.names, platform, repcond1, repcond2) {
  
  # Step 1: Initialization
  long <- length(condi.names)  # Number of conditions
  filt <- list()  # List to store filtered data frames
  
  # Step 2: Filtering and assignment
  if (platform == 1 | platform == 2){
    for (i in 1:long) {
      filt[[i]] <- raw[raw[, condi.names[i]] > 0, ]
    }
  } else if (platform == 3){
    for (i in 1:long) {
      filt[[i]] <- raw[complete.cases(raw[, condi.names[i]]), ]
    }
  }
  
  # Step 3: Counting proteins
  proteins <- sapply(filt, nrow)
  
  # Step 4: Visualization
  barplot(
    proteins,
    main = "Proteins identified",
    xlab = "Samples",
    ylab = "Total proteins",
    cex.names = 0.5,
    names.arg = condi.names,
    col = c(rep("light green", repcond1), rep("light blue", repcond2)),
    las = 2,
    ylim = c(0, max(proteins) + 200)
  )
  
}

#Venn Diagram

venn_diagram <- function(df, unique_proteins, label1, label2, color1, color2){
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
  
  venn.plot <- venn.diagram(
    x = x,
    category.names = c(label1, label2),
    filename = NULL,
    output = TRUE,
    imagetype = "png",
    main = "Venn Diagram",
    fill = c(color1, color2)
  )
  return(venn.plot)
}


##################################################################################
#Filtering

filter_valids <- function(df, unique_proteins, conditions1, min_count, at_least_one = FALSE, LOG2.names, labeltype) {
  
  
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
  
  if (labeltype == 1){
    
  
    cond.filter <- sapply(1:length(cond.names), function(i) { #En nuestro caso se itera del 1 al 2
      df2 <- common_df[cond.names[[i]]]   # Extrae las columnas de inter?s
      df2 <- as.matrix(df2)   # Lo convierte en matriz
      sums <- rowSums(is.finite(df2)) # Cuenta el numero de valores validos para cada condicion 
      sums >= min_count[i]   # Calculates whether min_count requirement is met si el n?mero de datos es igual o superior al elemento en min_count se establece true
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
    
    condition1_names <- cond.names[[1]] #Manejamos grupos 
    condition2_names <- cond.names[[2]]
    
    for (i in 1:nrow(common_df)) { #Probar a poner el CV multiplicado por 100 al resultado y no al denominador
      valuesA <- as.numeric(common_df[i, condition1_names])
      valuesB <- as.numeric(common_df[i, condition2_names])
      valuesC <- c(valuesA, valuesB)
      cv1 <- (sd(valuesA, na.rm = TRUE) / mean(valuesA, na.rm = TRUE)) * 100
      cv2 <- (sd(valuesB, na.rm = TRUE) / mean(valuesB, na.rm = TRUE)) * 100
      cv3 <- (sd(valuesC, na.rm = TRUE) / mean(valuesC, na.rm = TRUE)) * 100
      common_df$CV_Control[i] <- cv1
      common_df$CV_Tratamiento[i] <- cv2
      common_df$CV_TratamientoyControl[i] <- cv3
    
      }
  
  } else if (labeltype ==2){
    common_df <- df
    common_df <- remove_missing(common_df)
    
    condition1_names <- cond.names[[1]] #Manejamos grupos 
    condition2_names <- cond.names[[2]]
    
    for (i in 1:nrow(common_df)) { #Probar a poner el CV multiplicado por 100 al resultado y no al denominador
      valuesA <- as.numeric(common_df[i, condition1_names])
      valuesB <- as.numeric(common_df[i, condition2_names])
      valuesC <- c(valuesA, valuesB)
      cv1 <- (sd(valuesA, na.rm = TRUE) / mean(valuesA, na.rm = TRUE)) * 100
      cv2 <- (sd(valuesB, na.rm = TRUE) / mean(valuesB, na.rm = TRUE)) * 100
      cv3 <- (sd(valuesC, na.rm = TRUE) / mean(valuesC, na.rm = TRUE)) * 100
      common_df$CV_Control[i] <- cv1
      common_df$CV_Tratamiento[i] <- cv2
      common_df$CV_TratamientoyControl[i] <- cv3
      
    }
  }
  
  
  return(common_df)  
}

##################################################################################
#Normalization 

median_centering <- function(df, LOG2.names) {
  
  df[, LOG2.names] <- lapply(LOG2.names, #A las columnas de inter?s, se procede a calcular el LOG2-mediana de cada 
                             function(x) {
                               LOG2 <- df[[x]] #Se asigna el valor calculado
                               LOG2[!is.finite(LOG2)] <- NA   # Excluimos los missing values del c?lculo de la mediana 
                               gMedian <- median(LOG2, na.rm = TRUE)
                               LOG2 - gMedian #Se calcula el valor.
                             }
  )
  
  return(df)
}

normalization_func <- function(df, LOG2.names, method){
  if (method == "trimMean"){
    df[, LOG2.names] <- normalizeThis(dat = df[, LOG2.names], method = method, trimFa = 0.4)
  } else if (method != "trimMean"){
    df[, LOG2.names] <- normalizeThis(dat = df[, LOG2.names], method = method)
  }
  return(df)
}



##################################################################################
#Imputation

impute_data <- function(df, LOG2.names, width = 0.3, downshift = 1.8) {
  
  impute.names <- sub("LOG2", "impute", LOG2.names)
  
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
  
  impute.names <- sub("LOG2", "impute", LOG2.names)
  df[impute.names] <- lapply(LOG2.names, function(x) !is.finite(df[, x])) #Creamos columnas donde se ve si imputar o no
  
  #Imputaci?n
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
  y <- na.omit(y)
  A <- rowMeans(y, na.rm = TRUE) #Hago las medias de las filas
  CV <- (matrixStats::rowSds(data.matrix(y), na.rm = TRUE)/A)^2 #Calculo de la desviaci?n estandar relativa invocando el paquete matrixStats y usando su funci?n rowSds
  res <- data.frame(mean = A, CV = CV)
  plot(A, CV ,  ylim = c(min(CV)-0.001, max(CV) +0.001),  ...)
  if(trend){ 
    fit <- limma::loessFit(CV, A) #Invocamos la funcion loess del paquete limma, para hacer una regresi?n local
    o <- order(A)
    lines(A[o], fit$fitted[o], lwd =2, col = "red")
  }
  
  return(res)
}

boxplot_function <- function(df, cond.names,  ...){
  df <- na.omit(df[cond.names])
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

scatterplot_function <- function(df, colname1, colname2){
  corr_plot <- ggplot(df, aes(df[,colname1], df[,colname2])) +
    geom_point() +
    theme_minimal() +
    labs(y = colname2, x = colname1)
  
  ggMarginal(corr_plot, type = "densigram")
  
}

qqplot_function <- function(df, colname1, colname2, color){

  p <- data.frame(qqnorm(x = df[,colname1], y = df[,colname2]))
  ggplot(p, aes(x, y)) +
    geom_point(pch = 21,
               colour = color) +
    theme_minimal() +
    geom_qq_line(aes(sample = y)) +
    ylab("Sample dist") +
    xlab("Theoretical dist") + 
    ggtitle("Q-Q plot")
}

corrplot_function <- function(df, display, tl.col = "black", addCoef.col = "black"){
  df_cor <- cor(df, use = "complete.obs", method = "pearson") 
  corrplot(df_cor, method=display, tl.col = "black", addCoef.col = "black")
}

##################################################################################
#Differential expression analysis
##################################################################################
#Limma function
statistical_analysis <- function(df, LOG2.names,test, paired = FALSE, replicas_condicion1, replicas_condicion2, condition1, condition2, logfcup, logfcdown, sig, adjval, statval, unique_proteins, way){
  
  if (test == 2){
    condition1_names <- grep(condition1, LOG2.names, value = TRUE)
    condition2_names <- grep(condition2, LOG2.names, value = TRUE)
    
    #Control columns
    for (i in 1:replicas_condicion1){
      nam <- paste("control_sample", i, sep = "")
      assign(nam, condition1_names[i])
    }  #LOG2.WT1,#LOG2.WT2,#LOG2.WT3,#LOG2.WT4
    
    #Treatment columns
    for (i in 1:replicas_condicion2){
      nam <- paste("prob_column", i, sep = "")
      assign(nam, condition2_names[i])
    }  #LOG2.WT_H2O2_1, #LOG2.WT_H2O2_2, #LOG2.WT_H2O2_3, #LOG2.WT_H2O2_4
    ct <- c()
    for (i in ls()[grep("control_sample", ls())]){
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
    logFC <- fit.eb$coefficients[, 2] #Calculo del log fold-change
    p.value <- fit.eb$p.value[, 2]    # p-valor moderado correspondiente al estad?stico t moderado.
    adj.P.Val <- p.adjust(p.value, method = adjval)
    
    if (statval == 1){
      expression <- case_when(logFC >= logfcup & -log10(p.value) >= -log10(sig) ~ "Up-regulated",
                              logFC <= logfcdown & -log10(p.value) >= -log10(sig) ~ "Down-regulated",
                              TRUE ~ "Unchanged")#labels para expresion 
    } else if (statval == 2){
      expression <- case_when(logFC >= logfcup & -log10(adj.P.Val) >= -log10(sig) ~ "Up-regulated",
                              logFC <= logfcdown & -log10(adj.P.Val) >= -log10(sig) ~ "Down-regulated",
                              TRUE ~ "Unchanged")
    }
    results.eb <- data.frame(logFC, p.value, adj.P.Val, expression)
    results_rownames <-rownames(results.eb) #Obtenemos el nombre de las filas
    Protein <- c()                          #Creamos un vector vacio
    Protein_description <- c()
    for (i in results_rownames){
      new_value <- df[i, "Protein"]     #iteramos los rownames en la columna deseada
      Protein <- c(Protein, new_value)      #Creamos la columna con identificadores
      new_desc <- df[i, "Protein_description"]
      Protein_description <- c(Protein_description, new_desc) #Creamos la columna con descriptores. 
    }
    results.eb$Protein <- Protein 
    results.eb$Protein_description <- Protein_description
    row.names(results.eb) <- make.unique(results.eb$Protein)
    
    if (way == 1){
      return(results.eb)
    } else if (way == 2){
      #Empezamos a trabajar con las unicas control
      unique.control <- as.data.frame(unique_proteins[1])
      fusioncontrol.df <- data.frame(matrix(ncol = 6))
      n <-  c("logFC", "p.value", "adj.P.Val", "expression", "Protein", "Protein_description")
      colnames(fusioncontrol.df) <- n
      
      for (i in unique.control$Protein){
        new <- c(min(results.eb$logFC, na.rm = TRUE)-2, min(results.eb$p.value, na.rm = TRUE)/100, min(results.eb$adj.P.Val, na.rm = TRUE)/100, "Down-regulated", i, unique.control$Protein_description[unique.control$Protein==i] )
        fusioncontrol.df <- rbind(fusioncontrol.df, new)
      }
      fusioncontrol.df <- fusioncontrol.df[-1,]  
      
      #Empezamos a trabajar con la unicas tratamiento
      unique.treatment <- as.data.frame(unique_proteins[2])
      fusiontreatment.df <- data.frame(matrix(ncol = 6))
      n <-  c("logFC", "p.value", "adj.P.Val", "expression", "Protein", "Protein_description")
      colnames(fusiontreatment.df) <- n
      
      for (i in unique.treatment$Protein){
        new <- c(max(results.eb$logFC, na.rm = TRUE)+2, min(results.eb$p.value, na.rm = TRUE)/100, min(results.eb$adj.P.Val, na.rm = TRUE)/100, "Up-regulated", i, unique.treatment$Protein_description[unique.treatment$Protein==i] )
        fusiontreatment.df <- rbind(fusiontreatment.df, new)
      }
      fusiontreatment.df <- fusiontreatment.df[-1,] 
      
      
      unique.proteins.limma <- rbind(fusioncontrol.df, fusiontreatment.df)
      row.names(unique.proteins.limma) <- unique.proteins.limma$Protein
      
      
      unique.proteins.limma[c("logFC", "p.value", "adj.P.Val")] <- sapply(unique.proteins.limma[c("logFC", "p.value", "adj.P.Val")], as.numeric) 
      
      results.eb <- rbind(unique.proteins.limma, results.eb)
      row.names(results.eb) <- make.unique(results.eb$Protein)
      return(results.eb)
    }
    
  } else if (test == 1){
    condition1_names <- grep(condition1, LOG2.names, value = TRUE)
    condition2_names <- grep(condition2, LOG2.names, value = TRUE)
    # Assuming df is your data frame with duplicate row names
    
    for (i in 1:nrow(df)) {
      valuesA <- as.numeric(df[i, c(condition1_names)])
      valuesB <- as.numeric(df[i, c(condition2_names)])
      
      if (sum(is.finite(valuesA)) >= 2 && sum(is.finite(valuesB)) >= 2) {
        testResults <- t.test(x = valuesA, y = valuesB, paired = paired)
        
        df$pValue[i] <- testResults$p.value
        df$tStat[i] <- testResults$statistic
        df$logFC[i] <- mean(valuesB, na.rm = TRUE) - mean(valuesA, na.rm = TRUE)
        results.eb <- data.frame(df$logFC, df$pValue)
        
      } else {
        df$pValue[i] <- NA
        df$tStat[i] <- NA
        df$logFC[i] <- NA
      }
      
      
      
    }
    colnames(results.eb) <- c("logFC", "p.value")
    results.eb$adj.P.Val <- p.adjust(results.eb$p.value, method = adjval)
    if (statval==1){
      results.eb$expression <- case_when(results.eb$logFC >= logfcup & -log10(results.eb$p.value) >= -log10(sig) ~ "Up-regulated",
                                         results.eb$logFC <= logfcdown & -log10(results.eb$p.value) >= -log10(sig) ~ "Down-regulated",
                                         TRUE ~ "Unchanged")#labels para expresion 
    } else if (statval==2){
      results.eb$expression <- case_when(results.eb$logFC >= logfcup & -log10(results.eb$adj.P.Val) >= -log10(sig) ~ "Up-regulated",
                                         results.eb$logFC <= logfcdown & -log10(results.eb$adj.P.Val) >= -log10(sig) ~ "Down-regulated",
                                         TRUE ~ "Unchanged")
    }
    
    
    Protein <- c()                          #Creamos un vector vac?o
    Protein_description <- c()
    results_rownames <- rownames(results.eb)
    for (i in results_rownames){
      new_value <- df[i, "Protein"]     #iteramos los rownames en la columna deseada
      Protein <- c(Protein, new_value)      #Creamos la columna con identificadores
      new_desc <- df[i, "Protein_description"]
      Protein_description <- c(Protein_description, new_desc) #Creamos la columna con descriptores. 
    }
    results.eb$Protein <- Protein 
    results.eb$Protein_description <- Protein_description
    row.names(results.eb) <- make.unique(results.eb$Protein)
    
    if (way == 1){
      return(results.eb)
    } else if (way == 2){
      #Empezamos a trabajar con las unicas control
      unique.control <- as.data.frame(unique_proteins[1])
      fusioncontrol.df <- data.frame(matrix(ncol = 6))
      n <-  c("logFC", "p.value", "adj.P.Val", "expression", "Protein", "Protein_description")
      colnames(fusioncontrol.df) <- n
      
      for (i in unique.control$Protein){
        new <- c(min(results.eb$logFC, na.rm = TRUE)-2, min(results.eb$p.value, na.rm = TRUE)/100, min(results.eb$adj.P.Val, na.rm = TRUE)/100, "Down-regulated", i, unique.control$Protein_description[unique.control$Protein==i] )
        fusioncontrol.df <- rbind(fusioncontrol.df, new)
      }
      
      fusioncontrol.df <- fusioncontrol.df[-1,]  
      
      #Empezamos a trabajar con la unicas tratamiento
      unique.treatment <- as.data.frame(unique_proteins[2])
      fusiontreatment.df <- data.frame(matrix(ncol = 6))
      n <-  c("logFC", "p.value", "adj.P.Val", "expression", "Protein", "Protein_description")
      colnames(fusiontreatment.df) <- n
      
      for (i in unique.treatment$Protein){
        new <- c(max(results.eb$logFC, na.rm = TRUE)+2, min(results.eb$p.value, na.rm = TRUE)/100, min(results.eb$adj.P.Val, na.rm = TRUE)/100, "Up-regulated", i, unique.treatment$Protein_description[unique.treatment$Protein==i] )
        fusiontreatment.df <- rbind(fusiontreatment.df, new)
      }
      fusiontreatment.df <- fusiontreatment.df[-1,] 
      
      
      unique.proteins.limm <- rbind(fusioncontrol.df, fusiontreatment.df)
      row.names(unique.proteins.limm) <- unique.proteins.limm$Protein
      
      unique.proteins.limm[c("logFC", "p.value", "adj.P.Val")] <- sapply(unique.proteins.limm[c("logFC", "p.value", "adj.P.Val")], as.numeric) 
      
      results.eb <- rbind(unique.proteins.limm, results.eb)
      row.names(results.eb) <- make.unique(results.eb$Protein)
      return(results.eb)
    }
###    
    
  }
  
}


##################################################################################
#Plots

volcano_plot <- function(limma, title, label, statval){
  
  #top <- readline("Introduzca cuantas proteinas quiere etiquetar:")
  if (statval == 1){

    logFC <- limma$logFC
    p.value <- limma$p.value
    protein_ids <- limma$Protein
    
    
    rx <- c(-1, 1)*max(abs(logFC))*1.1 #Para definir nuestro eje X cogemos el m?ximo de la distribuci?n
    ry <- c(0, ceiling(max(-log10(p.value), -log10(p.value)))) #Para definir nuestro eje Y cogemos el el m?nimo y el m?ximo de la distribuci?n
    
    plot <- plot_ly(data = limma, x = ~logFC, y = ~-log10(p.value), text = protein_ids,
                    type = "scatter", mode = "markers",
                    color = ~expression, colors = c("green3", "gray78", "firebrick3"),
                    size = I(label))
    
    plot <- plot %>% layout(title = title,
                            xaxis = list(title = list(text ='Log2 Fold Change')),
                            yaxis = list(title = list(text = 'Log10 p-value')))
    
    return(plot)
  } else if (statval ==2){

    logFC <- limma$logFC
    adj.P.Val <- limma$adj.P.Val
    protein_ids <- limma$Protein
    
    
    rx <- c(-1, 1)*max(abs(logFC))*1.1 #Para definir nuestro eje X cogemos el m?ximo de la distribuci?n
    ry <- c(0, ceiling(max(-log10(limma$adj.P.Val), -log10(limma$adj.P.Val)))) #Para definir nuestro eje Y cogemos el el m?nimo y el m?ximo de la distribuci?n
    
    rx <- c(-1, 1)*max(abs(logFC))*1.1 #Para definir nuestro eje X cogemos el m?ximo de la distribuci?n
    ry <- c(0, ceiling(max(-log10(adj.P.Val), -log10(adj.P.Val)))) #Para definir nuestro eje Y cogemos el el m?nimo y el m?ximo de la distribuci?n
    
    plot <- plot_ly(data = limma, x = ~logFC, y = ~-log10(adj.P.Val), text = protein_ids,
                    type = "scatter", mode = "markers",
                    color = ~expression, colors = c("green3", "gray78", "firebrick3"),
                    size = I(label))
    
    plot <- plot %>% layout(title = title,
                            xaxis = list(title = list(text ='Log2 Fold Change')),
                            yaxis = list(title = list(text = 'Log10 p-value')))
    
    return(plot)
  }
  
}


pca <- function(x, group, rep1, rep2){
  x <- na.omit(x[group])
  tdata<-t(x[group])
  
  pc<-prcomp(tdata)
  eigenValues <- pc$sdev ^2
  totalEigenValue <- sum(eigenValues)
  pcContribution <- 100 * eigenValues / totalEigenValue
  contribution1 <- round(100 * eigenValues[1] / totalEigenValue)
  contribution2 <- round(100 * eigenValues[2] / totalEigenValue)
  my_colors <-  rep(c("salmon", "blue"), c(rep1, rep2))
  
  plot(pc$x[,1], pc$x[,2], pch=16, cex = 3,  col = my_colors, main="Principal Component Analysis", xlab=paste0("PC1 [", contribution1, "%]"), ylab=paste0("PC2 [", contribution2, "%]"))
  text(pc$x[,1], pc$x[,2],labels=rownames(pc$x),pos=3,offset=0.4,cex=0.7)
}



my_heatmap <- function(data, cond.names, title){
  
  heatmap <- heatmap(as.matrix(data[cond.names]), labRow = data$Protein, main = title, Rowv = NULL, Colv = NA, col =greenred(75), cexCol = 0.6)
  legend(x="topleft", legend=c("2", "0", "-2"),fill=c("red", "black", "green"), title = "Log2FC")

  
}

my_heatmap_differential <- function(limma, data, cond.names, title){
  top_proteins <- bind_rows(   #Nos quedamos con las prote?nas sobre e infra expresadas en las comparativas WT
    limma %>% 
      filter(limma$expression == 'Up-regulated'),  
    limma %>% 
      filter(limma$expression == 'Down-regulated') 
  )
  
  diferential <- data.frame()
  for (i in top_proteins$Protein){   #Nos quedamos con la entrada del dataframe del dataset original de cada proteina de inter?s
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

dotplot_func <- function(terms, ...){
  
  dotplot(terms[[1]],  ...) +  facet_grid(.~Conditions)
  
}

gostplot_func <- function(terms, ...){
  gostplot(terms[[3]], ...) 
  
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

##################################################################################
#Interaction network analysis

interactions_up <- function(df, taxonid, score){
  #Subset de up-regulated y down-regulated
  string_db <- STRINGdb$new(version="11.5", species=taxonid, score_threshold=score, input_directory="", protocol="http")

  
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

interactions_down <- function(df, taxonid, score){
  string_db <- STRINGdb$new(version="11.5", species=taxonid, score_threshold=score, input_directory="", protocol="http")
  
  down_regulated <- subset(df, df$expression == "Down-regulated")
  
  down_mapped <- string_db$map(down_regulated, "Protein", removeUnmappedRows = TRUE)
  
  par(mfrow=c(1,1))
  hits_down <- down_mapped$STRING_id[1:100] 
  string_db$plot_network(hits_down)
  
  
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