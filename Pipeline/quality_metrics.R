venn_diagram <- function(df, unique_proteins, ...){
  #Proteinas en muestras control
  cond1_set <- df%>%
    filter(df$Protein != unique_proteins[[2]]$Protein[1])
  for (i in unique_proteins[[2]]$Protein){
    cond1_set <- cond1_set%>%
      filter(cond1_set$Protein != i)
  }
  
  #Proteinas en muestras tratamiento
  cond2_set <- df%>%
    filter(df$Protein != unique_proteins[[1]]$Protein[1])
  for (i in unique_proteins[[1]]$Protein){
    cond2_set <- cond2_set %>%
      filter(cond2_set$Protein != i)
  }
  
  label1 <- readline("Introduzca nombre de condicion control:")
  label2 <- readline("Introduzca nombre de condicion tratamiento:")
  
  x <- list( cond1_set$Protein,  cond2_set$Protein)
  names(x) <- c(label1, label2)
  ggvenn(x, ...)
  ggsave("vennDiagram.tiff", units="in", width=7, height=7, dpi=300)
}

plotCV2 <- function(y, trend = TRUE, main= "Imputation check", ...){

  A <- rowMeans(y, na.rm = TRUE) #Hago las medias de las filas
  filename <- "Imputation check.tiff"
  tiff(filename = filename, units="in", width=5, height=4, res=300)
  CV <- (matrixStats::rowSds(data.matrix(y), na.rm = TRUE)/A)^2 #Calculo el CV invocando el paquete matrixStats y usando su función rowSds
  print(summary(CV))
  res <- data.frame(mean = A, CV = CV)
  plot(A, CV ,  ylim = c(min(CV)-0.001, max(CV) +0.001),  ...)
  if(trend){ 
    fit <- limma::loessFit(CV, A) #Invocamos la funcion loess del paquete limma, para hacer una regresión local
    o <- order(A)
    lines(A[o], fit$fitted[o], lwd =2, col = "red")
    dev.off()
  }
  
  return(res)
}

boxplot_function <- function(df, cond.names, ...){
  filename <- "Boxplot of normalized intensities.tiff"
  tiff(filename = filename, units="in", width=5, height=4, res=300)
  par(mfrow=c(1,1), font.lab=2, cex.lab=1, font.axis=2, cex.axis=1, cex.main=1, las = 1)
  boxplot(df[, cond.names], names = cond.names,  ylim = c(min(df[,cond.names])-1, max(df[,cond.names]) +1), main="Boxplot normalized Intensities", ylab = "Intensities", las=2, ...)
  dev.off()
}

histogram <- function(df, cond.names, color){
  dir.create("./histograms")
  histDirectory <- paste0(finalDirectory,"/histograms")
  setwd(histDirectory)
  
  for (i in cond.names){
    temp_name <- paste0("Histogram of normalized intensities", i)
    filename <- paste0(temp_name, ".tiff")
    tiff(filename = filename, units="in", width=5, height=4, res=300)
    hist(as.numeric(df[,i]), main = i, xlab = 'Intensity-value', col = color)
    dev.off()
  }
  
}

imputation_state <- function(df.F, df.FNI, cond.names){
  filename1 <- "Proportion of missing_1.tiff"
  tiff(filename = filename1, units="in", width=5, height=4, res=300)
  par(mfrow=c(2,2))
  aggr(df.F[, cond.names], delimiter = "NA", labels = names(df.F[cond.names]), cex.axis = .5)
  dev.off()
  filename2 <- "Proportion of missing_2.tiff"
  tiff(filename = filename2, units="in", width=5, height=4, res=300)
  aggr(df.FNI[, cond.names], delimiter = "NA", labels = names(df.F[cond.names]), cex.axis = .5)
  dev.off()
}