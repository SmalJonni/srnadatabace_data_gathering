rule combine_knowledges:
	input:
		lit="resources/Literature.tsv",
		db="results/Database_knowledge/combined.tsv",
		ngs="result/cmbn/denovo.tsv"
	output:
		"results/combined/aggregated.tsv"
	conda:
		"../envs/filter.yaml"
	threads:
		1
	notebook:
		"../notebooks/aggregate_all.r.ipynb"

rule remove_empty_sequences:
	input:
		rules.combine_knowledges.output
	output:
		"results/filter/nonempty_seq.tsv"
	conda:
		"../envs/filter.yaml"
	threads:
		1
	notebook:
		"../notebooks/rm_empty_seq.r.ipynb"
	
rule tax2taxid:
	input:
		"results/filter/nonempty_seq.tsv"
	output:
		"results/cleaning/taxmap.tsv"
	conda:
		"../envs/filter.yaml"
	threads:
		1
	notebook:
		"../notebooks/tax2id.r.ipynb"

rule integrate_cleans:
	input:
		pre_clean="results/filter/nonempty_seq.tsv",
		tax_cleandict="results/cleaning/taxmap.tsv"
	output:
		"results/cleaning/cleaned_cmbn.tsv"
	conda:
		"../envs/filter.yaml"
	threads:
		1
	notebook:
		"../notebooks/integrate.r.ipynb"

rule collapse_duplicates:
	input:
		"results/cleaning/cleaned_cmbn.tsv"
	output:
		"results/filter/cleaned.tsv"
	conda:
		"../envs/filter.yaml"
	threads:
		1
	notebook:
		"../notebooks/collapse_duplicates.r.ipynb"
