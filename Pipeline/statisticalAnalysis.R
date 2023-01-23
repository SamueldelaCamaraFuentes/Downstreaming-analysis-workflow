limma.analysis <- function(df, ...){
  orden <- readline("Introduzca como quiere ordenar su tabla resultante: teclee 1 para p-value o 2 para q-value:")
  paired <- readline("Introduzca 1 si su experimento no es pareado o 2 si lo es:")
  #Columnas control
  for (i in 1:replicas_condicion1){
    nam <- paste("control", i, sep = "")
    assign(nam, readline("Introduzca el nombre de su columna control:"))
  }  #LOG2.WT1,#LOG2.WT2,#LOG2.WT3,#LOG2.WT4
  
  #Columnas problema
  for (i in 1:replicas_condicion2){
    nam <- paste("prob_column", i, sep = "")
    assign(nam, readline("Introduzca el nombre de su columna tratamiento:"))
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
  
  if (paired == 1){
    control <- rep(1, replicas_condicion1)
    treatment <- rep(2, replicas_condicion2)
    design <- model.matrix(~factor(c(control, treatment)))
  } else if (paired == 2){
    pairinfo = factor(rep(1:replicas_condicion1,2))
    control <- rep(1, replicas_condicion1)
    treatment <- rep(2, replicas_condicion2)
    design <- model.matrix(~pairinfo+factor(c(control, treatment)))
  }
  log_FCup <- readline("Cota superior del logFC significativo:")
  log_FCdown <- as.numeric(readline("Cota inferior del logFC significativo:"))
  sig <- readline("Introduzca un valor significativo de p-value:")
  
  dat <- df[, c(ct, tr)]
  n <- dim(dat)[1]
  fit <- lmFit(dat, design)
  fit.eb <- eBayes(fit, ...)
  print(colnames(fit.eb))
  logFC <- fit.eb$coefficients[, 2] #Cálculo del log fold-change
  p.mod <- fit.eb$p.value[, 2]    # p-valor moderado correspondiente al estadístico t moderado.

  expression <- case_when(logFC >= log_FCup & -log10(p.mod) >= sig ~ "Up-regulated",
                          logFC <= log_FCdown & -log10(p.mod) >= sig ~ "Down-regulated",
                          TRUE ~ "Unchanged")#labels para expresion 
  adj.P.Val <- p.adjust(p.mod, method = "fdr")
  results.eb <- data.frame(logFC, p.mod, adj.P.Val, expression)
  if (orden == 1){
    results.eb <- results.eb[order(results.eb$p.mod), ] #Ordenamos por p.mod
  } else if (orden ==2){
    results.eb <- results.eb[order(results.eb$q.mod), ]
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
  write_xlsx(results.eb, paste0(finalDirectory,"/Analisis_expresion_diferencial.xlsx"))
  return(results.eb)
}