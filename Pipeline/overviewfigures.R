volcano_plot <- function(limma, conditions){
  
  top <- readline("Introduzca cuantas proteinas quiere etiquetar:")
  logFC <- limma$logFC
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
    scale_color_manual(values = c("green3", "gray78", "firebrick3")) +
    guides(colour = guide_legend(override.aes = list(size=1.5))) +
    ggtitle("Volcano plot",conditions) +
    labs(y = "-log10 p-value", x = "log2 Fold change") +
    geom_label(data = top_proteins,
               mapping = aes(logFC, -log10(p.mod), label = Protein),
               size = 1.5)
  ggsave("volcano.tiff", units="in", width=7, height=7, dpi=300)
  
  
}

pca <- function(x, group){
  filename <- "PCA_plot.tiff"
  tiff(filename = filename, units="in", width=9, height=9, res=300)
  kmc<-length(unique(group))
  tdata<-t(x[group])
  km <-kmeans(tdata,(kmc/4),2000) 
  pc<-prcomp(tdata)
  plot(pc$x[,1], pc$x[,2], pch=16, col = c("salmon", "blue"),main="K-means Clustering of PCA", xlab="PC1", ylab="PC2")
  text(pc$x[,1], pc$x[,2],labels=rownames(pc$x),pos=3,offset=0.4,cex=0.7)
  #pc<-cbind(pc$x[,1], pc$x[,2])
  #ordispider(pc, factor(km$cluster), label = TRUE)
  #ordihull(pc, factor(km$cluster), lty = "dotted")
  dev.off()
}


my_heatmap <- function(limma, data, cond.names, conditions, whole = FALSE, ...){
  #@limma: es la salida de la función limma.fit para la comparativa de replicas WT
  #@data: es el dataset filtrado, normalizado e imputado. 
  
  top_proteins <- bind_rows(   #Nos quedamos con las proteínas sobre e infra expresadas 
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
  
  
  #plots
  #length <- length(LOG2.names)
  filename = "Heatmap of differentially expressed proteins.tiff"
  tiff(filename = filename, units="in", width=14, height=9, res=300)
  heatmap <- heatmap(as.matrix(diferential[cond.names]),  labRow = diferential$Protein, main =conditions, ...)
  legend(x="bottomright", legend=c("2", "0", "-2"),fill=c("red", "black", "green"), title = "Log2FC")
  dev.off()
  
  if (whole == TRUE){
    filename2 = "Heatmap of proteins.tiff"
    tiff(filename = filename2, units="in", width=14, height=9, res=300)
    heatmap <- heatmap(as.matrix(data[cond.names]), labRow = data$Protein, main = "Heatmap completo", ...)
    legend(x="bottomright", legend=c("2", "0", "-2"),fill=c("red", "black", "green"), title = "Log2FC")
    
  }
  
  dev.off()
  return(top_proteins)
}