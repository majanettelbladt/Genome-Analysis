import pandas as pd
import matplotlib.pyplot as plt
import os

# Ange sökvägen till din mapp med HTSeq-count-filer
path = "/Users/majanettelbladt/Desktop/RNA-Seq"

files = [f for f in os.listdir(path) if f.endswith(".txt")]

all_counts = []

for file in files:
    file_path = os.path.join(path, file)

    df = pd.read_csv(file_path, sep="\t", header=None, names=["gene", "count"])
    df = df[~df["gene"].str.startswith("__")]
    df["count"] = pd.to_numeric(df["count"])

    all_counts.extend(df["count"].tolist())

# Create histogram
plt.figure(figsize=(10,6))
plt.hist(all_counts, bins=100, log=True, color="mediumseagreen", edgecolor="black")
plt.title("Histogram over read counts")
plt.xlabel("Read count per gene")
plt.ylabel("Number of genes (log-skala)")
plt.tight_layout()
plt.savefig("read_count_histogram.png")
plt.close()

# Statistik
all_counts_series = pd.Series(all_counts)
total_genes = len(all_counts_series)
expressed_1 = (all_counts_series > 0).sum()
expressed_10 = (all_counts_series >= 10).sum()
expressed_100 = (all_counts_series >= 100).sum()

print(f"Totalt number of genes: {total_genes}")
print(f"Genes with >0 counts: {expressed_1} ({expressed_1/total_genes:.2%})")
print(f"Genes with ≥10 counts: {expressed_10} ({expressed_10/total_genes:.2%})")
print(f"Genes with ≥100 counts: {expressed_100} ({expressed_100/total_genes:.2%})")
print("Histogram saved as: read_count_histogram.png")
