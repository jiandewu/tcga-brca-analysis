# TCGA Breast Cancer Analysis

This repository contains analysis code for TCGA breast cancer studies, with a
focus on triple-negative breast cancer (TNBC), RNA-seq, miRNA, and Illumina
HumanMethylation450 data.

It was extracted from the former `workcodes` monorepo so that the breast cancer
work can be maintained independently. Manuscript drafts and generated binary
artifacts are intentionally excluded from this public repository.

## Repository layout

```text
.
├── analysis/
│   ├── Brca_Methy/            # TCGA 450K methylation workflow
│   ├── FPKMDElimmaedge.R      # RNA-seq differential-expression analysis
│   ├── miRNADElimmaedge.R     # miRNA differential-expression analysis
│   └── *.shed                 # legacy analysis command records
├── data_acquisition/
│   └── brca.shed              # legacy GDC acquisition commands
├── legacy/
│   ├── Fokhrul/               # historical BLIA/BLIS analyses
│   ├── MachineLearn/          # historical machine-learning commands
│   ├── RNAseq_counts/         # read-count preparation and analysis
│   ├── methylation/           # earlier methylation workflows
│   └── miRNA/                 # earlier miRNA workflow
├── environment.yml
└── .gitignore
```

## Environment setup

Create the ARM64-compatible Conda environment:

```bash
conda env create -f environment.yml
conda activate tcga-brca-analysis
R -q -e 'install.packages("riskyr", repos="https://cloud.r-project.org")'
```

`riskyr` is installed from CRAN because it is not available as a Conda package
for the platform used to validate this repository.

## Historical workflow order

The newer TCGA methylation workflow was designed to run in numerical order:

```text
analysis/Brca_Methy/1DataDownload.R
analysis/Brca_Methy/2Preparation.R
analysis/Brca_Methy/3DMCs.R
analysis/Brca_Methy/4DMRs.R
analysis/Brca_Methy/5DMRs.R
```

Other files are retained as historical analysis code rather than a single
end-to-end workflow. Review inputs and parameters before running them.

## Data and publication policy

Raw GDC downloads, controlled or patient-level data, intermediate matrices,
generated plots, and manuscript drafts are not stored in Git. Keep controlled
data in institutionally approved storage. Public-data manifests, checksums,
source code, and small non-sensitive result summaries may be committed when
appropriate.

## Portability status

This is an initial history-preserving extraction. Many legacy scripts still
contain absolute paths under `/home/jiande/` and depend on reference files or
intermediate data that are not included. These paths must be replaced with
project-root configuration before the workflows are considered portable or
reproducible.

The `legacy/` directory intentionally preserves older approaches for reference.
Its contents should be reviewed and migrated selectively rather than treated as
a supported production pipeline.
