# AMM
`AMM` is a command-line tool to run the Abstract Mediation Model. For details on the method and application, please see Weiner et al 2021, biorxiv. Please cite that manuscript if you use this command-line tool in your own work.

## Overview

`AMM` estimates two quantities:
1) p(k): The fraction of heritability mediated by the kth-closest gene
2) e(A): Mediated heritability enrichment of a gene set

What you need to provide for estimation:
1) GWAS summary statistics
2) Gene sets for either 1) learning p(k) or 2) estimation of e(A)

## Getting started

You'll need to download this repository to access AMM.py, which contains the functions to run `AMM`. To download this repository, run:
```  
git clone https://github.com/danjweiner/AMM21.git
```
The repository contains the Conda environment file (amm.yml) for running `AMM`. Note that this environment includes Python 2.7, which is required by the `AMM` dependency [LD Score Regression](https://github.com/bulik/ldsc). On that note, you will need to download `LDSC` to run `AMM`, which you can do at the previous link.

## Overview of `AMM` workflow

![alt text](https://github.com/danjweiner/AMM21/blob/main/AMM_reference_files/AMM_overview.png)
