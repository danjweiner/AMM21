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

We've broken down the `AMM` workflow into the following steps, which we will walk through one by one. In broad strokes, the workflow involves:

**Module 1-3**: Preparing for analysis by creating an annotation matrix based on your gene set of interest

**Module 4, 6.2, 7**: Estimating p(k). You can skip this branch all together if you already have p(k) estimates you'd like to use. For example, you could use our p(k) estimates from the manuscript (see below). 

**Module 5, 6.1, 8**: Estimated e(A), mediated heritability enrichment.

![alt text](https://github.com/danjweiner/AMM21/blob/main/AMM_reference_files/AMM_overview.png)

With that, let's get into it:

## Module 1: Estimate kn-matrix

In this step, we need to define the k-th closest gene to each SNP. In `AMM`, we've written a function that takes in gene positions, SNP positions, and outputs a SNP x gene rank matrix for each chromosome, where the elements of the matrix are gene names. Since this is a bit memory intensive, and since most users of `AMM` won't need to re-make these, **we've published ours that you can download so you don't have to do Module 1 yourself** [here](https://drive.google.com/drive/u/1/folders/1t3tSCILksoGI16zqc-4vvfRcMfMNXLH1). We created these matrices out to the 50th closest gene.

If you want to create your own, run the following `AMM` command:



