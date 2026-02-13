# Pipeline de traitements des données GranzymeB

Ce pipeline traite des données RNA-seq brutes (fichiers FASTQ) pour produire des fichiers BAM alignés, filtrés et triés, ainsi que des statistiques de qualité et de mapping. Il utilise les environnements FastQC, Trimmomatic, STAR et Samtools, activés via Micromamba.


## Description des variables 

Le fichier config.json contient toutes les variables dont on a besoin dans ce worflow. 
En premier lieu nous définissons les chemin d'accès vers nos espaces de travail. 
Ensuite, nous appelons les fichiers Fasta des chromosomes de référence. 
Puis, viennent les noms des échantillons à analyser.
Suivis 
Enfin, nous définissons les paramètres qui seront utilisés lors du trimming pour pouvoir les modifier plus facilement.  


## Les règles du Workflow

Ces règles sont définies dans le fichier CutNTag_QC.smk. Pour chaque règle, est précisé dans le titre, l'environnement avec lequel nous travaillons.


### Etape 1 – FastQC sur les fastq bruts

```bash
rule Fastqc_raw
```
Objectif : Vérifier la qualité initiale des séquences brutes avant tout traitement.
Nous obtenons en sortie un fichier html pour R1 et R2 contenant des graphiques pour détecter les problèmes tels que des séquences de faible qualité, une taille anormale de reads, des adaptateurs parasites.


### Etape 2 – Trimming avec Trimmomatic

```bash
rule Trimmomatic
```
Objectif : Supprimer les séquences d’adaptateurs et les bases de faible qualité pour améliorer par la suite la qualité de l’alignement.
Nous utilisons les fonctions : 
- ILLUMINACLIP : supprime les adaptateurs connus
- LEADING/TRAILING : supprime les bases de faible qualité au début et à la fin
- MINLEN : supprime les reads trop courts après trimming
Nous obtenons alors à nouveau des fichiers fastq de meilleurs qualité pour l'alignement.


### Etape 3 - FastQC après trimming

```bash
rule Fastqc_trim
```
Objectif : Vérifier que le trimming a amélioré la qualité des séquences.
Nous obetenons en sortie des fichiers html de FastQC pour R1 et R2 après trimming. Ils permettent de nous assurer que les bases de mauvaise qualité ont été retirées et que les adaptateurs ont été supprimés.


### Etape 4 – Alignement avec Star

```bash
rule star
```
Objectif : Mapper les reads nettoyés sur le génome de référence (hG38).
Le fichier BAM que l'on obtient en sortie contient les positions des reads alignés sur le génome. En parallèle, le fichier log est un contrôle qualité de l'alignement. Il fournit les métriques de l'alignement (pourcentage de reads mappés, multimapping, etc).


### Etape 5 - Samtools stats

```bash
rule samtools_stats
```
Objectif : Générer des statistiques détaillées sur le fichier BAM produit par l'alignement.
samtools stats calcule des informations comme la couverture, le nombre de reads, les taux de duplication…
plot-bamstats produit des graphiques ensuite des graphiques regroupant ces informations.


### Etape 6 - Filter Star Bam

```bash
rule Filter_Bam
```
Objectif : Retirer les reads de faible qualité.
Seuls les reads avec un score de mapping suffisement élevés sont conservés. Cela permet d'avoir des données fiables pour les analyses.


### Etape 7 - Sort Star Bam

```bash
rule Sort_Bam
```
Objectif : Filtrer le fichier BAM pour les analyses ultérieures.
Le BAM est trié par coordonnées génomiques grâce à samtools sort. 
samtools index crée un index pour accéder rapidement aux régions spécifiques.


## Les modifications 

```bash
02/02/2026 Isaline git commit -m "suppression de toutes les variables dans config.json qu'on ne va pas utiliser"
05/02/2026 Isaline git commit -m "première modification et création des règles correspondant aux 3 premières étapes du pipeline (fast qc sur les fichiers bruts, trimming, fast qc sur les fichiers trimmés)"
```

