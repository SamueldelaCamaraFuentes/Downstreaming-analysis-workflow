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
    }
}

