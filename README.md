# advanced-python
This project performs a comprehensive analysis of gene expression data to identify differentially expressed genes (DEGs) between tumor and normal samples. The analysis involves loading and cleaning data from three files: "Gene_Expression_Data.xlsx", "Gene_Information.csv", and "Sample_Information.tsv". Sample names in the gene expression data are updated based on their phenotypes, followed by splitting the data into tumor and normal groups. The average gene expression for each group is calculated, and fold changes are determined to identify genes with significant expression differences. Genes with a fold change magnitude greater than 5 are flagged, and their upregulation or downregulation is noted. The exploratory data analysis (EDA) includes visualizations such as histograms, bar charts, heatmaps, and clustermaps to provide insights into gene expression patterns and distributions across chromosomes and sample types. Findings show significant differences in gene expression between tumor and normal samples, with notable upregulation in tumors.
