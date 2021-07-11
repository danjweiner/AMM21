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

`AMM` is written to be run on a cluster, which allows you to submit an array job to speed things up! You'll note that many modules require a flag `--iterator`, where you pass the task array environment variable.  

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

```
python amm.py\
	--m 1\
	--iterator [iterator variable; submit 1 job per autosome, so jobs: 1-22]\
	--genes_ref [text file with 1 row per gene, with a column of gene names and a column of gene locations in BP, by chromosome, ending .txt]\
	--snp_ref_bim [a text file mapping SNP RSID to SNP BP location; a .bim file from 1KG or equivalent works well, by chromosome, ending .bim]\
	--kn_k [How many rank columns you would like in the output matrix; we used 50]\
	--snp_loc_col [Column index of SNP location (0-indexed)]\
	--genes_loc_col [Column index of gene location (0-indexed)]\
	--out [Directory where you would like the kn-matrix to print out]\
	--snp_rsid_col [Column index of SNP name (0-indexed)]\
	--genes_ensg_col [Column index of Gene name (0-indexed)]
```

More concretely, the command might look like:

```
python path_to_amm/amm.py\
	--m 1\
	--iterator "$SGE_TASK_ID"\
	--genes_ref path_to_reference_files/gnomad_gene_location_guide_17661_chr\
	--snp_ref_bim path_to_reference_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.\
	--kn_k 50\
	--snp_loc_col 3\
	--genes_loc_col 7\
	--out path/kn_matrix/\
	--snp_rsid_col 1\
	--genes_ensg_col 1
```
The full file names are `path_to_reference_files/gnomad_gene_location_guide_17661_chr[#].txt` and `path_to_reference_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.[#].bim`

## Module 2: LD score annotations

The goal of this step is to create a per-chromosome SNP-by-proximity annotation matrix. The entries in this matrix are 1 or 0, denoting whether the kth closest gene to the SNP is in a gene set. Thus, for each gene set, you will create a new set of 22 annotation matrices, one per chromosome. 

You need two ingredients to make an annotation matrix: 

The first is the kn-matrix from Module 1, which you either created or downloaded premade. 

The second ingredient are gene sets. Each gene set is described by a list of genes in a text file, one gene per line, no header. Some important notes:
* Our premade kn-matrices use ENSG gene IDs (i.e. ENSG00000130164), so your text file with gene names must have ENSG gene names). 
* If you have a gene in your gene set that isn't in the kn-matrix, its as if it doesn't exist. So, its critical to cross-reference your gene set with the genes used to create the kn-matrix. If you use your kn-matrices, cross reference with `gnomad_gene_location_guide_17661.txt` in the AMM_reference_files folder in this repository.
* Each gene set gets its own text file, ending in .txt. Place these files in the AMM working directory 

Now you're ready to make annotations:

```
python amm.py\
	--m 2\
	--iterator [iterator variable; submit 1 job per autosome, so jobs: 1-22]\
	--set_names [gene set summary file, see below]\
	--kn_in_dir [directory where your kn-matrices live]\
	--out [AMM working directory]\
	--kn_k [number of proximity annotations you would like; must be less than or equal to the number of columns in the kn-matrix]
```

More concretely, the command might look like:

```
python path_to_amm/amm.py\
	--m 2\
	--iterator "$SGE_TASK_ID"\
	--set_names amm_gs.txt\
	--kn_in_dir path/kn_matrix/\
	--out amm_working_directory/\
	--kn_k 50
```
A few additional notes:
* `amm_gs.txt` lists your gene set files, one per line. The gene set files are listed *without* extensions and *without* .txt. For instance, if in your AMM working directory you have a `gene_set_1.txt` and `gene_set_2.txt`, then your file amm_gs.txt would be:

```
gene_set_1
gene_set_2
```
## Module 3: 
