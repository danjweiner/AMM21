from __future__ import division
from __future__ import print_function
import numpy as np
import pandas as pd
import sys
import os
import argparse 
import os.path

#Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("--snp_ref_bim", help="a text file mapping SNP RSID to SNP BP location; a .bim file from 1KG or equivalent works well, by chromosome, ending .bim", required=False)
parser.add_argument("--snp_loc_col", help="Column index of SNP location (0-indexed)", required=False, type = int)
parser.add_argument("--snp_rsid_col", help="Column index of SNP name (0-indexed)", required=False, type = int)
parser.add_argument("--genes_ref", help="Text file with 1 row per gene, with a column of gene names and a column of gene locations in BP, by chromosome, ending .txt", required=False)
parser.add_argument("--genes_loc_col", help="Column index of gene location (0-indexed)", required=False, type = int)
parser.add_argument("--genes_ensg_col", help="Column index of Gene name (0-indexed)", required=False, type = int)
parser.add_argument("--kn_k", help="Module 1: Number of rank columns you would like in the output matrix; Module 2: number of columns in your kn matrix", required=False, type = int)
parser.add_argument("--kn_k_out", help="Module 2: number of proximity annotations you would like; must be less than or equal to the number of columns in the kn-matrix", required=False, type = int)
parser.add_argument('--m', help="module number", required=True, type = int)
parser.add_argument('--out', help="AMM working directory", required=True)
parser.add_argument('--kn_in_dir', help="directory where your kn-matrices live", required=False)
parser.add_argument('--set_names', help="Path with file containing set names, one per line, no quotes, copied to working directory, ending must be .txt but excluded from this file", required=False)
parser.add_argument('--lds_ref_binary', help="path to LD reference panel, .bim/.bed/.fam format, up to the chr integer", required=False)
parser.add_argument('--lds_snps_out', help="Path to list of SNPs to restrict LD score output to, ends .snps", required=False)
parser.add_argument('--pk_vector', help="p(k) list, same ordering as the --pk_size bunch list, one entry per bin, a list of floats (no comma separation)", required=False, type=float, nargs='+')
parser.add_argument('--pk_size', help="list specifying bunching; one integer per group, number of integers is number of groups. Sum of integers in list must equal number of proximity annotations so all annotations end up in a group, (no comma separation)", required=False, type=int, nargs='+')
parser.add_argument('--ss_list', help="List of summary statistics names and paths for analysis, in format [trait_name]:[trait path]", required=False)
parser.add_argument('--control_list', help="Control files for analysis, in format [control_name]:[control path], where control path ends before the chr integer", required=False)
parser.add_argument('--freq_file', help="Allele frequencies, available from LDSC repository, ending before chr integer", required=False)
parser.add_argument('--weights', help="LD score regression weights, available from LDSC repository, ending before chr integer", required=False)
parser.add_argument('--ldsc_path', help="Path to ld score regression", required=False)
parser.add_argument('--control_name', help="Name of control annotation in the input files", required=False)
parser.add_argument('--control_name_m', help="Path to control annotation in the reference .M file, up to chromosome integer", required=False)
parser.add_argument('--pk_out_name', help="Name of file created in p(k) meta-analysis", required=False)
parser.add_argument("--n_genes_total", help="Number of genes in genome", required=False, type = int)
parser.add_argument('--n_jk', help="Number of jackknife blocks in LD score regression, defaults to 200", nargs='?', const=200, type=int, default=200, required = False)
parser.add_argument("--iterator", help="Iterator for module; see Github Readme for details", required=False, type = int)
parser.add_argument("--which_regression", help="Which LD score regression mode: 1 for imputed pk regression for enrichment, 2 for bunched regression to allow p(k) estimation", required=False, type = int)
args = parser.parse_args()

#Functions
def check_ascending(file, location_column):
    return(file.sort_values(location_column, ascending = True).reset_index(drop=True))
    
def difference_matrix(snps, genes, snp_location_column, gene_location_column):
    n_columns = len(genes)
    n_index = len(snps)
    
    range_columns = np.arange(n_columns)
    range_index = np.arange(n_index)
    
    gene_location = pd.Series(genes[gene_location_column])
    snp_location = pd.Series(snps[snp_location_column])
    
    diff_mtx = abs(pd.DataFrame(np.array([snp_location]*n_columns).T, index=range_index, columns=range_columns) 
               - pd.DataFrame(np.array([gene_location]*n_index), index=range_index, columns=range_columns))
    
    return(diff_mtx)
                              
def k_closest(difference_mtx, k_genes, genes, ensg_column, snps, rsid_col):
    
    k_closest = {}

    for i in range(len(difference_mtx)):
        a = difference_mtx.iloc[[i,]].T
        a['gene'] = a.index
        a = a.sort_values(by=[i])
        a = a.reset_index(drop=True)
        a = a['gene'].head(k_genes)
        k_closest[i] = a
    
    k_closest_mtx = pd.DataFrame(k_closest).T
    k_closest_mtx_ensg = k_closest_mtx.replace(pd.Series(genes[ensg_column], index=genes.index).to_dict())
    k_closest_mtx_ensg.index = snps[rsid_col]
    
    return(k_closest_mtx_ensg)
    
def filter_k_closest(k_closest_mtx, ensg_gene_list):
	return(k_closest_mtx.isin(ensg_gene_list).astype(int))

def reduce_k(matrix, k):
	return(matrix.iloc[:,0:k])

def lds_from_pk(pk_array, group_length_array, lds_file):
    return(pd.concat([lds_file.iloc[:, 0:3], 
           pd.DataFrame((lds_file.iloc[:, 3:]*np.repeat(pk_array, group_length_array)).sum(axis = 1))], axis=1))

def m_5_50_median(m_5_50_file):
    return(pd.DataFrame(m_5_50_file.median(axis = 1).round()))

def bunch_annot(annotation_file):
    return((pd.DataFrame(annotation_file.sum(axis = 1)) > 0).astype(int))

def jk_enrichment(tau0, tauA, total_genes, set_genes):
    return((tau0 + tauA)/(tau0 + ((set_genes/total_genes)*tauA)))
    
def enrichment_parameters(table, n_jk, trait_name, gs_name, set_size, total_size):
    df = pd.DataFrame([table.mean(), table.std()*np.sqrt(n_jk-1)]).T
    df.columns = ['Enrichment_mean', 'Enrichment_se']
    df['Enrichment_z'] = (df['Enrichment_mean']-1)/df['Enrichment_se']
    df = df.sort_values(by = 'Enrichment_z',ascending = False)
    df['Trait'] = trait_name
    df['Gene_Set'] = gs_name
    df['Gene_set_size'] = set_size
    df['Total_size'] = total_size
    return(df)

def lds_for_pk(lds_array, bunch_array):
    sum_results = {}
    counter = 1
    left_bound = 3
    for bunch in bunch_array:
        right_bound = left_bound + bunch
        sum_results[counter] = pd.DataFrame(lds_array.iloc[:, left_bound:right_bound].sum(axis = 1))
        left_bound = right_bound
        counter = counter + 1
    return(pd.concat([lds_array.iloc[:, 0:3], pd.concat(sum_results, axis=1).sum(axis=1, level=0)], axis=1))

def m_5_50_median_to_pk(m_5_50_file, bunch_array):
    median_results = {}
    counter = 1
    left_bound = 0
    for bunch in bunch_array:
        right_bound = left_bound + bunch
        median_results[counter] = pd.DataFrame(m_5_50_file.iloc[:, left_bound:right_bound].median(axis = 1).round())
        left_bound = right_bound
        counter = counter + 1
        
    return(pd.concat(median_results, axis=1).sum(axis=1, level=0))

def bunch_annot_to_pk(annot_file, bunch_array):
    sum_max1_results = {}
    label = {}
    counter = 1
    left_bound = 0
    for bunch in bunch_array:
        right_bound = left_bound + bunch
        sum_max1_results[counter] = (pd.DataFrame(annot_file.iloc[:, left_bound:right_bound].sum(axis = 1)) > 0).astype(int)
        left_bound = right_bound
        counter = counter + 1
    return(pd.concat(sum_max1_results, axis=1).sum(axis=1, level=0))

def tauA_pk(coefficient_table, bunches_array):
    return(pd.DataFrame(coefficient_table.iloc[:,0:len(bunches_array)].dot(bunches_array)))

def pk(tauA_vector, tauK_vectors, jk_blocks):
    p_k_mean = {}
    p_k_sd = {}
    counter = 1
    for k in range(len(tauK_vectors)):
        pk_jk = tauK_vectors[counter]/tauA_vector
        p_k_mean[counter] = pk_jk.mean()
        p_k_sd[counter] = np.std(pk_jk)*np.sqrt(jk_blocks-1)
        counter = counter + 1
        
    results = pd.DataFrame([pd.concat(p_k_mean), pd.concat(p_k_sd)]).T.reset_index(drop=True)
    results.columns = ['p(k)', 'p(k)_se']
    
    return(results)

if __name__ == '__main__':

	print('Abstract Mediation Model')
	print('Command-line tool (c) Daniel Weiner, 2021')
	print('Version 0.1')
	print('')

	if args.m == 1: ##Module 1: KN matrix (SA: 1-22)
		print('Running module 1: K-nearest matrix')
		print('')
		
		snps_file = pd.read_csv(args.snp_ref_bim + str(args.iterator) + '.bim', header = None, sep = "\t")
		snps_file = check_ascending(snps_file, args.snp_loc_col)
		genes_file = pd.read_csv(args.genes_ref + str(args.iterator) + '.txt', header = None, skiprows=1, sep = " ")
		genes_file = check_ascending(genes_file, args.genes_loc_col)    

		k_closest(difference_matrix(snps_file, genes_file, args.snp_loc_col, args.genes_loc_col), args.kn_k, genes_file, args.genes_ensg_col, snps_file, args.snp_rsid_col).to_csv(args.out + "kn_matrix_k" + str(args.kn_k) + "_chr" + str(args.iterator) + '.kn.gz', sep = " ", header = False, index = True, compression='gzip')

		print('Completed module 1 for chromosome ' + str(args.iterator))

	if args.m == 2:	##Module 2: Annotation matrix (SA: 1-22)
		print('Running module 2: K-nearest matrix -> annotation matrix')
		print('')
		
		sets_names = pd.read_csv(args.set_names, header = None)[0].tolist()
		kn = pd.read_csv(args.kn_in_dir + "kn_matrix_k" + str(args.kn_k) + "_chr" + str(args.iterator) + '.kn.gz',  index_col = 0, header = None, sep = " ")
		del kn.index.name

		for index in range(len(sets_names)): 
			gene_set = pd.read_csv(args.out+sets_names[index]+".txt", header = None)[0].tolist()
			filter_k_closest(reduce_k(kn, args.kn_k_out), gene_set).to_csv(args.out+sets_names[index]+'_chr'+str(args.iterator)+'.annot.gz', sep = " ", header = True, index = False, compression='gzip')

		print('Completed module 2 for chromosome ' + str(args.iterator))
		
	if args.m == 3: #LD score (SA: gene sets * chromosomes)
		print('Running module 3: annotation matrix -> LD scores')
		print('')
		
		sets_names = pd.read_csv(args.set_names, header = None)[0].tolist()
		output = {}
		counter = 1
		for set in sets_names:
			for chrom in range(1,23):
				output[counter] = [set, chrom]
				counter = counter + 1
        
		lds_parameters = pd.DataFrame(output).T
		set_name = lds_parameters.iloc[int(args.iterator)-1,0]
		chrom = lds_parameters.iloc[int(args.iterator)-1,1]

		regression_command = args.ldsc_path + " --l2 --bfile " + args.lds_ref_binary + str(chrom) + " --ld-wind-cm 1 --annot " + args.out + set_name + "_chr" + str(chrom) +  ".annot.gz --thin-annot --out " + args.out + set_name + "_chr" + str(chrom) + " --print-snps " + args.lds_snps_out + str(chrom) + ".snps"
		os.system(regression_command)

		print('Completed module 3 for iteration ' + str(args.iterator))

	if args.m == 4: ##Module 4: Bunching of LD scores to allow ultimate calculation of p(k) (SA: gene sets * 22 chromosomes)
		print('Running module 4: LD scores -> bunched LD scores to allow ultimate calculation of p(k)')		
		print('')
		sets_names = pd.read_csv(args.set_names, header = None)[0].tolist()

		output = {}
		counter = 1
		for set in sets_names:
 		   for chrom in range(1,23):
 		       output[counter] = [set, chrom]
 		       counter = counter + 1
  		      
		lds_parameters = pd.DataFrame(output).T
		set_name = lds_parameters.iloc[int(args.iterator)-1,0]
		chrom = lds_parameters.iloc[int(args.iterator)-1,1]

		lds = pd.read_csv(args.out + set_name + "_chr" + str(chrom) + ".l2.ldscore.gz", sep = "\t")
		m550 = pd.read_csv(args.out + set_name + "_chr" + str(chrom) + ".l2.M_5_50", sep = "\t", header = None)
		annot = pd.read_csv(args.out + set_name + "_chr" + str(chrom) + ".annot.gz", sep = " ")

		lds_for_pk(lds, args.pk_size).to_csv(args.out + set_name + "_calculate_pk_chr" + str(chrom) + ".l2.ldscore.gz", compression = "gzip", sep = " ", index = None)
		m_5_50_median_to_pk(m550, args.pk_size).to_csv(args.out + set_name + "_calculate_pk_chr" + str(chrom) + ".l2.M_5_50", index = None, header = None, sep = " ")
		bunch_annot_to_pk(annot, args.pk_size).to_csv(args.out + set_name + "_calculate_pk_chr" + str(chrom) + ".annot.gz", index = None, sep = " ", compression='gzip')

		print('Completed module 4 for iteration ' + str(args.iterator))

	if args.m == 5: #Bunched LD scores for enrichment (and .M_5_50 and annot files) (SA: gene sets * 22 chromosomes)
		print('Running module 5: LD scores -> bunched LD scores')
		print('p(k) vector is: ' + str(args.pk_vector))
		print('bin vector is: ' + str(args.pk_size))
		print('')
		
		sets_names = pd.read_csv(args.set_names, header = None)[0].tolist()
		output = {}
		counter = 1
		for set in sets_names:
		    for chrom in range(1,23):
		        output[counter] = [set, chrom]
		        counter = counter + 1
        
		lds_parameters = pd.DataFrame(output).T
		set_name = lds_parameters.iloc[int(args.iterator)-1,0]
		chrom = lds_parameters.iloc[int(args.iterator)-1,1]

		lds = pd.read_csv(args.out + set_name + "_chr" + str(chrom) + ".l2.ldscore.gz", sep = "\t")
		m550 = pd.read_csv(args.out + set_name + "_chr" + str(chrom) + ".l2.M_5_50", sep = "\t", header = None)
		annot = pd.read_csv(args.out + set_name + "_chr" + str(chrom) + ".annot.gz", sep = " ")

		lds_from_pk(args.pk_vector, args.pk_size, lds).to_csv(args.out + set_name + "_pk_impute_" +  "chr" + str(chrom) + ".l2.ldscore.gz", compression = "gzip", sep = "\t", index = None)
		m_5_50_median(m550).to_csv(args.out + set_name + "_pk_impute_" +  "chr" + str(chrom) + ".l2.M_5_50", index = None, header = None)
		bunch_annot(annot).to_csv(args.out + set_name + "_pk_impute_" +  "chr" + str(chrom) + ".annot.gz", index = None, compression = "gzip")

		print('Completed module 5 for iteration ' + str(args.iterator))

	if args.m == 6: ##Module 6: Regression (SA: gene sets * traits * control models) 
		print('Running module 6: LD scores -> LD score regression (Note: Assumed LD score jackknife = 200 blocks, unless specified with optional flag)')
		print('')
		sets_names = pd.read_csv(args.set_names, header = None)[0].tolist()
		ss_list = pd.read_csv(args.ss_list, header = None)[0].tolist()
		control_list = pd.read_csv(args.control_list, header = None)[0].tolist()

		counter = 1
		reg_combo = {}

		for g in sets_names:
			for k in control_list:
				for n in ss_list:
					reg_combo[counter] = [g,k,n]
					counter = counter + 1

		results = pd.DataFrame(reg_combo).T
		results[['control_name', 'control_path']] = results[1].str.split(":",expand=True)
		results[['ss_name', 'ss_path']] = results[2].str.split(":",expand=True)
		results = results.drop(columns=[1,2])
		results.columns = ['set_name', 'control_name', 'control_path', 'ss_name', 'ss_path']

		set_name = results.iloc[int(args.iterator)-1,0]
		control_name = results.iloc[int(args.iterator)-1,1]
		control_path = results.iloc[int(args.iterator)-1,2]
		ss_name = results.iloc[int(args.iterator)-1,3]
		ss_path = results.iloc[int(args.iterator)-1,4] 

		if args.which_regression == 1:
			regression_command = args.ldsc_path + " --h2" + " " + ss_path + " --ref-ld-chr " + args.out + set_name + "_pk_impute_chr," + control_path + " --out " +  args.out + set_name + ".imputed_pk." + control_name + "." + ss_name + " --print-coefficients " + " --print-delete-vals " + " --w-ld-chr " + args.weights + " --frqfile-chr " + args.freq_file + " --overlap-annot " 
			os.system(regression_command)
			print('Completed module 6.1 (imputed p(k) for iteration ' + str(args.iterator))

		if args.which_regression == 2:
			regression_command = args.ldsc_path + " --h2" + " " + ss_path + " --ref-ld-chr " + args.out + set_name + "_calculate_pk_chr," + control_path + " --out " +  args.out + set_name + "." + "calc_pk" + "." + control_name + "." + ss_name + " --print-coefficients " + " --print-delete-vals " + " --w-ld-chr " + args.weights + " --frqfile-chr " + args.freq_file + " --overlap-annot " 
			os.system(regression_command)
			print('Completed module 6.2 (bunched to estimate p(k) for iteration ' + str(args.iterator))

	if args.m == 7: ##Module 7: Estimation of p(k)
		print('Running module 7: Bunched coefficients -> estimation of p(k) (Note: Assumed LD score jackknife = 200 blocks, unless specified with optional flag)')	
		print('')
		sets_names = pd.read_csv(args.set_names, header = None)[0].tolist()
		
		ss_list = pd.read_csv(args.ss_list, header = None)[0].tolist()
		ss_list_enrichment = pd.DataFrame(ss_list)
		ss_list_enrichment[['ss_name', 'ss_path']] = ss_list_enrichment[0].str.split(":",expand=True)
		
		results_tauA = {}
		counter = 1
		for sets in sets_names:
			for trait in ss_list_enrichment['ss_name']:
				coefficient_matrix = pd.read_csv(args.out+sets + ".calc_pk." + args.control_name + "." + trait + ".part_delete", header = None, sep = " ")
				results_tauA[counter] = tauA_pk(coefficient_matrix, args.pk_size)
				counter = counter + 1
		tau_A_sum = pd.DataFrame(pd.concat(results_tauA, axis = 1).sum(axis = 1))

		k_counter = 1
		k_results = {}
		for k in range(len(args.pk_size)):
			results_tauk = {}
			counter = 1
			for sets in sets_names:
				for trait in ss_list_enrichment['ss_name']:
					coefficient_matrix = pd.read_csv(args.out+sets + ".calc_pk." + args.control_name + "." + trait + ".part_delete", header = None, sep = " ")
					results_tauk[counter] = coefficient_matrix.iloc[:,k]
					counter = counter + 1
	
			k_results[k_counter] = pd.DataFrame(pd.concat(results_tauk, axis = 1).sum(axis = 1))
			k_counter = k_counter + 1

		pk(tau_A_sum, k_results, args.n_jk).to_csv(args.out + args.pk_out_name, index = False, sep = " ")
		
		print('Completed module 7: p(k)')
		
	if args.m == 8: ##Module 8: Enrichment of traits in gene sets 
		print('Running module 8: LD score regression -> AMM enrichment (Note: Assumed LD score jackknife = 200 blocks, unless specified with optional flag)')
		print('')
		
		m_files = {}
		for chrom in range(1, 23):
			m_files[chrom] = pd.read_csv(args.control_name_m+"."+str(chrom)+".l2.M", sep = "\t", header = None)
		m_files_sum_across_chr = pd.concat(m_files).sum()
		
		sets_names = pd.read_csv(args.set_names, header = None)[0].tolist()
		
		ss_list = pd.read_csv(args.ss_list, header = None)[0].tolist()
		ss_list_enrichment = pd.DataFrame(ss_list)
		ss_list_enrichment[['ss_name', 'ss_path']] = ss_list_enrichment[0].str.split(":",expand=True)
		
		for set_name in sets_names:
		    set_n = len(pd.read_csv(args.out + set_name + ".txt", header = None))
		    results = {}
		    for trait in ss_list_enrichment['ss_name']:
		        check_file = os.path.isfile(args.out + set_name + ".imputed_pk." + args.control_name + "." + trait + ".part_delete")
			if check_file == False:
				continue
			if check_file == True:
				tauA = pd.read_csv(args.out + set_name + ".imputed_pk." + args.control_name + "." + trait + ".part_delete", header = None, sep = " ")
		        	tauA_array = tauA[0]
		        
		        	tau0 = tauA.iloc[:,1:]
        			tau0.columns = range(tau0.shape[1])
        			tau0_array = (tau0.dot(m_files_sum_across_chr))/m_files_sum_across_chr[0]
		
		        	results[trait] = enrichment_parameters(jk_enrichment(tau0_array, tauA_array, args.n_genes_total, set_n), args.n_jk, trait, set_name, set_n, args.n_genes_total)
    
		    	pd.concat(results, axis=0).sort_values(by = 'Enrichment_z', ascending = False).to_csv(args.out + set_name + "." + args.control_name + ".enrichment", header = True, index = None, sep = " ")
    
		print('Completed module 8: AMM enrichment')
