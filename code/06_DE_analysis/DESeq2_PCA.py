import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os

# Mapp där dina HTSeq-count filer finns
htseq_dir = '/Users/majanettelbladt/Desktop/RNA-Seq'  # ändra till din sökväg

# Läs in alla HTSeq-count filer (*.txt eller *.counts)
files = sorted(glob.glob(os.path.join(htseq_dir, '*.txt')))

# Skapa en tom dict för att lagra counts
counts_dict = {}

for f in files:
    sample_name = os.path.basename(f).replace('.txt', '')  # exempelvis ERR1797969_coverage.txt -> ERR1797969_coverage
    df = pd.read_csv(f, sep='\t', header=None, names=['gene', sample_name], index_col=0)
    counts_dict[sample_name] = df[sample_name]

# Slå ihop alla count-kolumner till en DataFrame
counts_df = pd.concat(counts_dict.values(), axis=1)
counts_df.columns = counts_dict.keys()

# Ta bort rader som inte är gener (HTSeq lägger ofta in __no_feature osv)
counts_df = counts_df[~counts_df.index.str.startswith('__')]

# Filter: behåll gener med minst 10 reads totalt
counts_df = counts_df[counts_df.sum(axis=1) >= 10]

# Log2-transformering med pseudocount 1
log_counts = np.log2(counts_df + 1)

# PCA
pca = PCA(n_components=2)
pcs = pca.fit_transform(log_counts.T)

# Skapa DataFrame för PCA-resultat
pca_df = pd.DataFrame(pcs, columns=['PC1', 'PC2'], index=log_counts.columns)

# Metadata - här kan du skapa en DataFrame manuellt för dina prover och deras grupper
metadata = pd.DataFrame({
    'sample': pca_df.index,
    'condition': ['Serum', 'Serum', 'Serum', 'BH', 'BH', 'BH'],  # exempel, justera efter dina prover
    'replicate': ['rep1', 'rep2', 'rep3', 'rep1', 'rep2', 'rep3']  # exempel
})
metadata = metadata.set_index('sample')

# Slå ihop PCA med metadata
pca_df = pca_df.join(metadata)

# Procent varians som varje komponent förklarar
explained_var = pca.explained_variance_ratio_ * 100

# Plot
plt.figure(figsize=(8,6))
sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='condition', style='replicate', s=100)
plt.xlabel(f'PC1 ({explained_var[0]:.1f}% variance)')
plt.ylabel(f'PC2 ({explained_var[1]:.1f}% variance)')
plt.title('PCA plot of RNA-seq samples based on HTSeq counts')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()
