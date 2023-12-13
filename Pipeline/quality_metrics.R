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

histogram <- function(df, cond.names, color){
  dir.create("./histograms")
  histDirectory <- paste0(finalDirectory,"/histograms")
  setwd(histDirectory)
  
  for (i in cond.names){
    temp_name <- paste0("Histogram of normalized intensities", i)
    filename <- paste0(temp_name, ".png")
   png(filename = filename, units="px", width=1000, height=1000, res=NA)
    hist(as.numeric(df[,i]), main = i, xlab = 'Intensity-value', col = color)
    dev.off()
  }
  
}

imputation_state <- function(df.F, df.FNI, cond.names, cex){
  filename1 <- "Proportion of missing_1.png"
  png(filename = filename1, units="px", width=1000, height=1000, res=NA)
  par(mfrow=c(2,2))
  aggr(df.F[, cond.names], delimiter = "NA", labels = names(df.F[cond.names]), cex.axis = .2)
  dev.off()
  filename2 <- "Proportion of missing_2.png"
  png(filename = filename2, units="px", width=1000, height=1000, res=NA)
  aggr(df.FNI[, cond.names], delimiter = "NA", labels = names(df.F[cond.names]), cex.axis = cex)
  dev.off()
}

scatterplot_function <- function(df, colname1, colname2){
  corr_plot <- ggplot(df, aes(df[,colname1], df[,colname2])) +
    geom_point() +
    theme_minimal() +
    labs(y = colname2, x = colname1)

  return(corr_plot)
  
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