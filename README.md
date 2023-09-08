
 ![Alt text](interface.gif) 
# Label Free Downstreaming Analysis 
En este repositorio se encuentra todo el código generado que comprende la herramienta de análisis downstreaming de datos de proteómica cuantitativa *Label free*.
<br>
<br>
En el presente repositorio podemos encontrar el software de análisis de datos de proteómica cuantitativa *Label free* enfocado en el análisis de datos tanto de DDA como de DIA procedentes de las plataformas computacionales MaxQuant, MSFragger y DIA-NN. Dentro de sus funcionalidades destacan el llevar a cabo un preprocesado de los datos, análisis de expresión diferencial, análisis funcional y de interacción de proteínas, junto con la visualización de cada una de las anteriores etapas mencionadas
Además en el presente repositorio se encuentra un *pipeline* para el análisis de datos de proteómica cuantitativa *label free* procedentes de MaxQuant cuyos archivos podemos encontrar en la carpeta "Pipeline" que contiene el archivo main.R eje conductor de las funcionalidades desarrolladas y los archivos dataprocessing.R quality_metrics.R, statisticalAnalysis.R, overviewfigures.R, enrichment.R y Network.R que contienen el conjunto de funciones implementadas. El conjunto de archivos presentes en la carpeta "Pipeline" desempeñan las funciones de la aplicación en RStudio.   
<br>
# Requisitos
- Archivo "ProteinGroups.txt" de MaxQuant. Encontramos un dataset de prueba en la carpeta "data".
- Archivo "combined_protein.tsv" de MSFragger.
- Archivo "report.pg.matrix.tsv" de DIA-NN.

<br>
# Outputs
1. Quality metrics:
- Boxplots
- Dispersion plot
- Scatter plots
- Pre/post imputation plots
- Histograms
- Principal Component Analysis
- Correlation plot
- Q-Q plot
2. Overview figures:
- Volcano plot
- Heatmap
3. Enrichment:
- Dotplot
- Barplot
- Manhattan plot
4. Interactions:
- Interaction networks
5. Data:
- *data set* tras preprocesamiento
- *data set* tras análisis estadístico
- *data set* tras análisis funcional
# Créditos
Paquetes que lo hacen posible: dplyr, ggplot2, ggvenn, VIM, gplots, gprofiler2, igraph, plotly, limma, clusterprofiler, enrichplot, stringdb.


