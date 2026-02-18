#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler

# Charger le fichier featureCounts
df = pd.read_csv("count_matrix.txt", sep="\t", comment="#")

# Colonnes de comptage
count_columns = df.columns[6:]
counts = df[count_columns]
counts.index = df["Geneid"]

# Supprimer gènes peu exprimés
counts = counts[counts.sum(axis=1) > 10]

# Garder les 1000 gènes les plus variables
top_var = counts.var(axis=1).sort_values(ascending=False).head(1000).index # On est pas obligées de garder les 1000 gènes les plus important --> à choisir 
counts = counts.loc[top_var]

# Normalisation Z-score par gène
scaler = StandardScaler()
counts_scaled = pd.DataFrame(
    scaler.fit_transform(counts.T).T,
    index=counts.index,
    columns=counts.columns
)

# Heatmap
plt.figure(figsize=(12, 50))
sns.heatmap(
    counts_scaled,
    cmap="vlag",
    center=0,
    xticklabels=True,
    yticklabels=True
)

plt.title("Gene-level Heatmap (Top 1000 variable genes)")
plt.tight_layout()
plt.savefig("gene_heatmap_1000.png", dpi=300)
plt.close()

print("Heatmap sauvegardée : gene_heatmap_1000.png")
