
 ![Alt text](interface.gif) 
# Label Free Downstreaming Analysis

In this repository, the downstreaming analysis code for quantitative proteomics analysis tool can be found.
<br>
<br>
In the present repository, we can find the quantitative proteomic data analysis software focused on analyzing both DDA and DIA label free and labelled data from MaxQuant, MSFragger, Proteome Discoverer and DIA-NN computational platforms. Its functionalities include preprocessing of data, differential expression analysis, functional and protein interaction analysis, along with visualization of each of the aforementioned stages.

Additionally, in this repository, there is a pipeline for analyzing label-free and labelled quantitative proteomic data from several computational platforms such as MaxQuant, MSFragger, Proteome Discoverer and DIA-NN, whose files can be found in the "Pipeline" folder containing the main.R file, which drives the developed functionalities, and the files dataprocessing.R, quality_metrics.R, statisticalAnalysis.R, overviewfigures.R, enrichment.R, and Network.R, which contain the set of implemented functions. The set of files in the "Pipeline" folder serves the functions of the application in RStudio.

<br>

# Requirement
- "ProteinGroups.txt" file from MaxQuant. a test dataset can be found on "data" folder.
- "combined_protein.tsv" file from MSFragger.
- "report.pg.matrix.tsv" file from DIA-NN.
- Output file from Proteome Discoverer. 

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
- *data set* after preprocessing
- *data set* after statistical analysis
- *data set* after functional analsysis.
  
# Credits
Packages that make it work: dplyr, ggplot2, ggvenn, VIM, gplots, gprofiler2, igraph, plotly, limma, clusterprofiler, enrichplot, stringdb.


