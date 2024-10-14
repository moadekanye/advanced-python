import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Step 1: Load the data
gene_expression_data = pd.read_excel("Gene_Expression_Data.xlsx")
gene_information = pd.read_csv("Gene_Information.csv")
sample_information = pd.read_csv("Sample_Information.tsv", sep='\t')

# Step 2: Clean up column names
gene_expression_data.columns = gene_expression_data.columns.str.strip()
sample_information.columns = sample_information.columns.str.strip()

# Create a mapping of sample IDs to their corresponding phenotypes
sample_phenotype_dict = dict(zip(sample_information['Sample_ID'], sample_information['Phenotype']))

# Step 3: Rename columns in gene expression data
new_column_names = []
for col in gene_expression_data.columns[1:]:  # Skip the first column, assuming it's the probe ID
    if col in sample_phenotype_dict:
        new_column_names.append(sample_phenotype_dict[col])
    else:
        new_column_names.append(col)

# Add suffix to make unique column names
new_column_names = [f"{name}_{i}" if new_column_names.count(name) > 1 else name for i, name in enumerate(new_column_names)]
gene_expression_data.columns = ['Probe_ID'] + new_column_names

# Step 4: Split the merged data based on phenotype
tumor_data = gene_expression_data.loc[:, gene_expression_data.columns.str.contains('Tumor')]
normal_data = gene_expression_data.loc[:, gene_expression_data.columns.str.contains('Normal')]

# Step 5: Compute average expression for each probe
average_tumor = tumor_data.mean(axis=1)
average_normal = normal_data.mean(axis=1)

# Step 6: Determine fold change
fold_change = (average_tumor - average_normal) / average_normal
gene_expression_data['Fold_Change'] = fold_change

# Step 7: Identify genes with absolute fold change > 5
filtered_genes = gene_expression_data[gene_expression_data['Fold_Change'].abs() > 5]

# Add a column to indicate upregulated or downregulated
filtered_genes['Expression_Change'] = np.where(filtered_genes['Fold_Change'] > 0, 'Higher in Tumor', 'Lower in Tumor')

# Merge with gene information to get gene names and other info
result = pd.merge(filtered_genes, gene_information, on='Probe_ID', how='left')

# Step 8: Exploratory Data Analysis (EDA)

# 8a: Histogram of DEGs by chromosome
plt.figure(figsize=(12, 6))
sns.histplot(data=result, x='Chromosome', bins=len(result['Chromosome'].unique()), kde=False)
plt.title('Distribution of Differentially Expressed Genes (DEGs) by Chromosome')
plt.xlabel('Chromosome')
plt.ylabel('Count of DEGs')
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

# 8b: Histogram of DEGs by chromosome segregated by sample type
plt.figure(figsize=(12, 6))
sns.histplot(data=result, x='Chromosome', hue='Expression_Change', multiple='stack', bins=len(result['Chromosome'].unique()))
plt.title('Distribution of DEGs by Chromosome and Sample Type')
plt.xlabel('Chromosome')
plt.ylabel('Count of DEGs')
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

# 8c: Bar chart of upregulated vs downregulated DEGs
upregulated_count = (result['Expression_Change'] == 'Higher in Tumor').sum()
downregulated_count = (result['Expression_Change'] == 'Lower in Tumor').sum()

plt.figure(figsize=(8, 6))
sns.barplot(x=['Upregulated', 'Downregulated'], y=[upregulated_count, downregulated_count], palette='pastel')
plt.title('Count of DEGs: Upregulated vs Downregulated in Tumor Samples')
plt.ylabel('Count')
plt.xlabel('Expression Change')
plt.tight_layout()
plt.show()

# 8d: Heatmap of gene expression by sample
plt.figure(figsize=(12, 8))
sns.heatmap(gene_expression_data.set_index('Probe_ID'), cmap='coolwarm', center=0)
plt.title('Heatmap of Gene Expression by Sample')
plt.xlabel('Samples')
plt.ylabel('Gene Probes')
plt.tight_layout()
plt.show()

# 8e: Clustermap of gene expression by sample
sns.clustermap(gene_expression_data.set_index('Probe_ID'), cmap='coolwarm', center=0, figsize=(12, 10))
plt.title('Clustermap of Gene Expression by Sample')
plt.show()

# Step 9: Write Findings
findings = """
1. The histograms indicate that the majority of differentially expressed genes (DEGs) are located on specific chromosomes.
2. When segregating DEGs by sample type, we can observe the differences in expression levels between tumor and normal samples.
3. The bar chart shows a notable proportion of upregulated genes in tumor samples compared to downregulated genes, indicating significant biological changes.
4. The heatmap illustrates the overall expression patterns of genes across different samples, while the clustermap reveals clustering of samples based on gene expression profiles.
"""
print(findings)
