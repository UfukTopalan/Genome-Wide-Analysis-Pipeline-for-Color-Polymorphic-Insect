# Genome-Wide Analysis Pipeline for a Color-Polymorphic Insect

Analysis pipeline, scripts, intermediate data, and results supporting:

> **Genome-Wide Analysis Reveals Altitude-Associated Divergence in a Color-Polymorphic Insect**
> Topalan & Sağlam (2026)

The full step-by-step pipeline (ANGSD/PCAngsd/NGSadmix/realSFS commands and R scripts, in
execution order) is documented in [`docs/Analyses_Pipeline_Topalan_Saglam_2026.pdf`](docs/Analyses_Pipeline_Topalan_Saglam_2026.pdf).

## Repository structure

```
.
├── docs/
│   └── Analyses_Pipeline_Topalan_Saglam_2026.pdf   # Full written pipeline (all steps, in order)
│
├── scripts/                # All R scripts used for plotting and downstream analysis
│   ├── plotPCA.R
│   ├── plotQQ.R
│   ├── run_Evanno_deltaK.R
│   ├── plotAdmix.R
│   ├── plotSFS.R
│   ├── plotFst.R
│   ├── run_mantel_test.R
│   ├── plotIBD.R
│   ├── plotIBDmulti.R
│   ├── plotDiversity.R
│   ├── plotMafTrend.R 
│   ├── cline_fitting_hzar_parboot.R
│   ├── summarize_clines.R
│   ├── plotDensity.R
│   └── plotHeatMapGeno.R
│
├── data/
│   ├── reference/           # Reference contigs (isophya_contigs_CAYMY.fasta) and BAM list
│   ├── sites/                # ANGSD -sites filter files (all-site, neutral, adaptive subsets)
│   ├── genotypes/            # Called genotypes (isophya71.geno.tsv)
│   └── metadata/             # Sample/population info (.info, .clst) and pairwise pop combinations
│
└── results/
    ├── admixture/            # NGSadmix likelihoods across K runs
    ├── fst/                  # Pairwise Fst tables (all sites, neutral, adaptive)
    ├── diversity/            # Theta/diversity statistics (all, neutral, adaptive)
    ├── association/          # Altitude-association MAF trend tables
    ├── cline/                # HZAR cline-fitting summary
    └── dapc/                 # DAPC discriminant scores, assignment probabilities, top SNPs
```

## Pipeline overview

See the PDF in `docs/` for exact commands, but in brief the workflow is:

1. **Master BEAGLE file generation** (ANGSD genotype likelihoods)
2. **PCA** (PCAngsd covariance matrix, eigen decomposition, SNP-loading QQ plots)
3. **Admixture** (NGSadmix, K = 1–5, Evanno's ΔK)
4. **Pairwise Fst / IBD** (folded SFS via realSFS, Mantel test)
5. **Partitioning neutral vs. adaptive loci** (HWE outliers, pcadapt)
6. **Neutral vs. adaptive diversity/Fst comparisons**
7. **Altitude association test** (ANGSD `-doAsso 2`)
8. **Allele-frequency trend plots** for altitude-associated loci
9. **Cline parameter estimation** (HZAR, parametric bootstrap)
10. **Genetic discrimination between colour morphs** (DAPC, Manhattan plot)

## Notes

- Large data files (`.fasta`, `.sites`, genotype `.tsv`) are tracked with **Git LFS** — see
  `.gitattributes`. Make sure Git LFS is installed before cloning/pulling
  (`git lfs install`), otherwise these files will appear as small pointer files.
- Sample/population codes throughout (e.g. `450`, `850`, ... `2300`) refer to collection
  site elevation in meters.
