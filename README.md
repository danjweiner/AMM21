# AMM
`AMM` is a command-line tool to run the Abstract Mediation Model. For details on the method and application, please see Weiner et al. 2021, bioRxiv. Please cite that manuscript if you use this command-line tool in your own work.

## Overview

`AMM` estimates two quantities:
1) p(k): The fraction of heritability mediated by the kth-closest gene
2) e(A): Mediated heritability enrichment of a gene set

The two main ingredients you'll need for estimation are:
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

In this step, we need to define the k-th closest gene to each SNP. In `AMM`, we've written a function that takes in gene positions, SNP positions, and outputs a SNP x gene rank matrix for each chromosome, where the elements of the matrix are gene names. Since this is a bit memory intensive, and since most users of `AMM` won't need to re-make these, **we've uploaded ours [here](https://drive.google.com/drive/u/1/folders/1t3tSCILksoGI16zqc-4vvfRcMfMNXLH1) so you don't have to run Module 1** . We created these matrices out to the 50th closest gene.

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
	--genes_ref /path_to_reference_files/gnomad_gene_location_guide_17661_chr\
	--snp_ref_bim /path_to_reference_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.\
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
	--set_names /path_to_set_names/amm_gs.txt\
	--kn_in_dir /path/kn_matrix/\
	--out /path_to_amm_working_directory/\
	--kn_k 50
```
A few additional notes:
* `amm_gs.txt` lists your gene set files, one per line. The gene set files are listed *without* extensions and *without* .txt. For instance, if in your AMM working directory you have a `gene_set_1.txt` and `gene_set_2.txt`, then your file amm_gs.txt would be:

```
gene_set_1
gene_set_2
```

## Module 3: Estimate LD scores

In this module, we use `LDSC` to estimate LD scores for the annotations we created in Module 2. The `AMM` call here makes it a bit easier to estimate LD scores for all of our gene sets at once. Refer to the `LDSC` documentation for more detail if needed. 

The setup is:

```
python amm.py\
	--m 3\
	--iterator [iterator variable; submit 1 job per autosome per gene set. For instance, if you want LD scores for 3 gene sets, you would submit 3*22 = 66 jobs]\
	--set_names [gene set summary file, see Module 2 for details]\
	--ldsc_path [path to `LDSC`]\
	--lds_ref_binary [path to LD reference panel, .bim/.bed/.fam format]\
	--lds_snps_out [path to list of SNPs per chromosome for which we will print LD scores]\
	--out [AMM working directory]\
```

More concretely, the command might look like:

```
python path_to_amm/amm.py\
	--m 3\
	--iterator "$SGE_TASK_ID"\
	--set_names /path_to_set_names/amm_gs.txt\
	--ldsc_path /path_to_ldsc/ldsc.py\
	--lds_ref_binary /path_to_reference_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.\
	--lds_snps_out /path_to_snp_files/snps.\
	--out /path_to_amm_working_directory/\
```
where the full file names are: `path_to_reference_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.[#].bim/bed/fam` and `path_to_snp_files/snps.[#].snps`. `amm_gs.txt` is the same file from Module 2.

**If you have pre-generated p(k) values and want to skip to estimating mediated heritability enrichment of a gene set, scroll down to Module 5** Pre-generated p(k) trained on constrained genes is available in the Supplementary Tables of the manuscript.

## Module 4: Bunch LD scores for estimating p(k)

The purpose of bunching LD scores is to estimate a mean p(k) for a group of annotations. Bunching means summing LD scores in the annotations being bunched into a group. For example, if in Module 3 you generated LD scores for the closest - 10th closest genes, perhaps instead of estimating p(k) for each rank (likely noisy), you can estimate the average for the closest, 2-5th, and 6th-10th closest genes. Module 4 allows you to specify which bins you would like to create:

```
python amm.py\
	--m 4\
	--iterator [iterator variable; submit 1 job per autosome per gene set. For instance, if you want LD scores for 3 gene sets, you would submit 3*22 = 66 jobs]\
	--out [AMM working directory]\
	--pk_size [list specifying bunching; one integer per group, number of integers is number of groups. Sum of integers in list must equal number of proximity annotations so all annotations end up in a group]\
	--set_names [gene set summary file, see Module 2 for details]
```
More concretely, the command might look like:
```
python path_to_amm/amm.py\
	--m 4\
	--iterator "$SGE_TASK_ID"\
	--out /path_to_amm_working_directory/\
	--pk_size 1 1 3 5 10 10 10 10\
	--set_names /path_to_set_names/amm_gs.txt
```
The annotations are bunched from start to end of the array. For example, the array given `1 1 3 5 10 10 10 10` sums to 50, which is the number of proximity annotations we created (i.e., closest - 50th closest genes across 50 columns, left to right). This list will create the following 8 bunches: closest gene, 2nd closest gene, 3rd-5th closest gene, 6th-10th closest gene, 11th-20th closest gene, 21st-30th closest gene, 31st-40th closest gene, 41st-50th closest gene.

## Module 6.2: LD Score Regression for estimating p(k)

Module 6 generates a regression call in `LDSC`, and submodule 6.2 is expecting the specific output generated in Module 4. This command will output regression coefficients for each of our annotations, which we will use next in Module 7 to estimate p(k). 

Module 6.2 is designed to allow you to run lots of regressions at once. Specifically, you provide (a) a list of at least 1 gene set (b) a list of at least 1 GWAS summary statistic (c) a list of at least 1 set of control annotations (like the Baseline LD model). Module 6.2 will run a `LDSC` regression for each combination; for example, if you provided 3 gene sets, 10 traits, and 1 control annotation, it would run 3 * 10 * 1 = 30 regressions. The outline of a Module 6.2 run is:

```
python amm.py\
	--m 6\
	--which_regression [specify which submodule to run]\
	--iterator [iterator variable; submit 1 job per regression. If you provided 3 gene sets, 10 traits, and 1 control annotation, submit 3 * 10 * 1 = 30 jobs]\
	--set_names [gene set summary file, see Module 2 for details]\
	--ss_list [list of trait names and summary statistics paths, see below for formating details]\
	--control_list [list of control annotations names and paths, see below for formating details]\
	--weights [LDSC regression weights]\
	--freq_file [LDSC allele frequency file]\
	--ldsc_path [path to LDSC]\
	--out [AMM working directory]
```
More concretely, the command might look like:
```
python path_to_amm/amm.py\
	--m 6\
	--which_regression 2\
	--iterator "$SGE_TASK_ID"\
	--set_names /path_to_set_names/amm_gs.txt\
	--ss_list /path_to_summary_statistics_file/amm_ss_full_47.txt\
	--control_list /path_to_control_list/control_list.txt\
	--weights /path_to_weights/weights.hm3_noMHC.\
	--freq_file /path_to_frequency_file/1000G.EUR.QC.\
	--ldsc_path /path_to_ldsc/ldsc.py\
	--out /path_to_amm_working_directory/
```
A few file format notes:

`path_to_summary_statistics_file/amm_ss_full_47.txt` is a text file listing the name and path to GWAS summary statistics files. Each line of the file requires format: [trait name]:[path to summary statistics file]. For example, you might have on one line: `standing_height:/path_to_summary_statistics/standing_height.txt`. The summary statistics format must be compatable with `LDSC`; see their documentation on [munge_sumstats.py](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation) for more details.

`/path_to_control_list/control_list.txt` is a text file listing the name and path to any control annotation you want to add to your regression. For example, it is common to add the [Baseline LD model](https://alkesgroup.broadinstitute.org/LDSCORE/) to `LDSC`. If you wanted to add this control annotation, you would in the format [control name]:[path to control LD scores], for example: `baseline_ld:/path_to_baseline_ld/baselineLD.` which will reference the files /path_to_baseline_ld/baselineLD.1.l2.ldscore.gz, /path_to_baseline_ld/baselineLD.2.l2.ldscore.gz, etc. 

The --weights and --freq_file arguments are standard from `LDSC` and can be downloaded from their repository.

# Module 7: Estimate p(k)

The purpose of Module 7 is to take the regression output from Module 6.2 and estimate p(k) from a single trait-gene set pair or (more commonly) through meta-analysis of more than 1 trait-gene set pair. The outline of a Module 7 run is:

```
python amm.py\
	--m 7\
	--out [AMM working directory]\
	--pk_size [list specifying annotation bunching size; same as Module 4]\
	--pk_control_name [name of control annotation in the input files]\
	--pk_out_name [outname of p(k) results file]\
	--set_names [gene sets you want to meta-analyze over]\
	--ss_list [list of summary statistics for meta-analysis, same format as in Module 6.2]
```
More concretely, the command might look like:
```
python path_to_amm/amm.py\
	--m 7\
	--out /path_to_amm_working_directory/\
	--pk_size 1 1 3 5 10 10 10 10\
	--pk_control_name baseline_ld\
	--pk_out_name cortex_liver_baseline_all_traits.txt\
	--set_names cortex_liver.txt\
	--ss_list /path_to_summary_statistics_file/amm_ss_full_47.txt
```
A few notes:

`--pk_control_name` We've written meta-analysis to be performed across only one set of control annotations at a time. If you created runs with different control annotations in Module 6.2, specify the name of the control annotation you want to include here. You select it by matching the argument of `--pk_control_name` with one of the names in `--control_list` from Module 6.2.

`--set_names` This is the same format as `--set_names` from Module 2. *You can specify a subset of gene sets to meta-analyze over* by editing this text file to include only a subset of gene set names.

`--ss_list` This is the same format as the summary statistics files above. In a similar way to `--set_names`, *you can specify a subset of traits to meta-analyze over* by subsetting the summary statistics text file to the traits you want.














# Software authorship
Daniel Weiner (Broad Institute)
