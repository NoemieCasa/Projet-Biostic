# FastQC

Ce repository contient les scripts pour télécharger, installer et executer les scripts nécéssaires au FastQC.

## Description des variables 


## Les règles du Workflow

### Etape 1 – Installer FastQC

```bash
sbatch scripts/Obtain_FastQC
```

## Les modifications 

```bash
02/02/2026 Isaline git commit -m "suppression de toutes les variables dans config.json qu'on ne va pas utiliser"
05/02/2026 Isaline git commit -m "première modification et création des règles correspondant aux 3 premières étapes du pipeline (fast qc sur les fichiers bruts, trimming, fast qc sur les fichiers trimmés)"
```

