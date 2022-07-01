from os import listdir

#  _____                             _                    
# |_   _|  __ _   _ _   __ _   ___  | |_                  
#   | |   / _` | | '_| / _` | / -_) |  _|                 
#   |_|   \__,_| |_|   \__, | \___|  \__|                 
#  ___                 |___/_        _     _              
# | _ \  _ _   ___   __| | (_)  __  | |_  (_)  ___   _ _  
# |  _/ | '_| / -_) / _` | | | / _| |  _| | | / _ \ | ' \ 
# |_|   |_|   \___| \__,_| |_| \__|  \__| |_| \___/ |_||_|
#                                       
#https://bioinf.shenwei.me/taxonkit/usage/#taxonkit                  
rule prepare_pytaxonkit_download:
	output:
		temp("results/setup/taxonkit/taxdump.tar.gz")
	shell:
		"wget -O {output} ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"

rule prepare_pytaxonkit_place:
	input:
		rules.prepare_pytaxonkit_download.output
	output:
		"results/setup/taxonkit/names.dmp",
		"results/setup/taxonkit/nodes.dmp",
		"results/setup/taxonkit/delnodes.dmp",
		"results/setup/taxonkit/merged.dmp"
	threads:
		1
	shell:
		"tar -zxf {input} -C $(dirname {output[0]} )"

checkpoint rule split_by_taxa:
	input:
		"results/filter/cleaned.tsv"
	output:
		directory("results/tax_split_fasta/")
	threads:
		1
	notebook:
		"../notebooks/tbl2fasta.py.ipynb"

def taxgatherer():
	k = listdir(checkpoints.split_by_taxa.get().output[0])
	return [i for i in k if i[0]!="."]


#Rule required due to inconvenience in sRNAFTarget
rule pickle_softlink:
	input:
		"resources/sRNARFTarget/PickledModelData/RFModel/sRNARFTargetModel.pickle"
	output:
		directory("PickledModelData")
	shell:
		"ln -s $(dirname $(dirname {input})) {output}"

checkpoint rule pull_mrnas:
	input:
		"results/filter/cleaned.tsv",
		rules.prepare_pytaxonkit_place.output
	output:
		directory("results/pulled_mrnas/")
	conda:
		"../envs/add_info.yaml"
	threads:
		1
	params:
		address=config["email_address"]
	notebook:
		"../notebooks/pull_mrnas.py.ipynb"

#Handle passing with wrapper is not possible in concurrent scenario due to output file being hardcoded
rule compute_targets:
	input:
		rules.pickle_softlink.output,
		m=lambda wildcard: checkpoints.pull_mrnas.get().output[0]+ wildcard+ "_mrna.fa",
		s="results/tax_split_fasta/{wc}_srna.fa"
#		m="results/pulled_mrnas/{wc}_mrna.fa"
	output:
		"results/targets/{wc}.tsv"
	params:
		pipeline="resources/sRNARFTarget/sRNARFTarget.nf"
	conda:
		"../envs/sRNARFTarget.yaml"
#	handover: 
#		True
	shadow:
		"shallow"
	threads:
		1
	shell:
		"""
		nextflow run resources/sRNARFTarget/sRNARFTarget.nf --s {input.s} --m {input.m} -w .
		mv  sRNARFTargetResult/Prediction_probabilities.csv  {output}
		"""


rule combine_targets:
	input:
		lambda wildcards:expand("results/targets/{wc}.tsv",wc=taxgatherer()),
		lambda wildcards:expand("results/tax_split_fasta/{wc}_srna.fa",wc=taxgatherer())
	conda:
		"../envs/add_info.yaml"
	output:
		"results/cmbnd/targets.tsv"
	conda:
		"../envs/add_info.yaml"
	threads:
		1
	notebook:
		"../notebooks/combine_targets.py.ipynb"

#  ___                              _                       
# / __|  ___   __   ___   _ _    __| |  __ _   _ _   _  _   
# \__ \ / -_) / _| / _ \ | ' \  / _` | / _` | | '_| | || |  
# |___/ \___| \__| \___/ |_||_| \__,_| \__,_| |_|    \_, |  
# / __| | |_   _ _   _  _   __  | |_   _  _   _ _   _|__/   
# \__ \ |  _| | '_| | || | / _| |  _| | || | | '_| / -_)    
# |___/  \__| |_|    \_,_| \__|  \__|  \_,_| |_|   \___|    
#
rule sequence_splitter:
        input:
                "results/filter/cleaned.tsv"
        output:
                "results/filter/sequences.tsv",
                "results/filter/main.tsv"
        conda:
                "../envs/add_info.yaml"
        threads:
                1
        notebook:
                "../notebooks/split.py.ipynb"
                                                           

rule tranform2fasta:
	input:
                "results/filter/sequences.tsv"
	output:
		"results/fasta/sequences.fa"
	conda:
		"../envs/add_info.yaml"
	notebook:
		"../notebooks/generate_fasta.py.ipynb"

rule viennaRNA:
	input:
		rules.tranform2fasta.output
	output:
		sec="results/viennarna/output_file.fs",
		ps=directory("results/viennarna/ps/")
	threads:
		128
	conda:
		"../envs/add_info.yaml"
	shadow: 
		"shallow"
	shell:
		"""
		mkdir -p  $(dirname {output.ps})
		RNAfold --infile={input} --jobs={threads} -t 4 --outfile=$(basename {output.sec})
		mkdir -p {output.ps}
		printf '%s\\0' *.ps | xargs -0 mv -t {output.ps}
		mv $(basename {output.sec} ) {output.sec}
		"""

def gatherer():
	k = listdir(checkpoints.split_hack.get().output[0])
	return [i for i in k if i[0]=="0"]

rule vienna_cmbn:
	input:
		rules.viennaRNA.output
	output:
		"results/cmbnd/vienna.tsv"
	threads:
		1
	conda:
		"../envs/add_info.yaml"
	notebook:
                "../notebooks/vienna_cmbn.py.ipynb"

rule lineage:
	input:
		"results/filter/main.tsv",
		rules.prepare_pytaxonkit_place.output
	output:
		"results/extened/lineage.tsv"
	threads:
		1
	conda:
		"../envs/add_info.yaml"
	notebook:
		"../notebooks/lineage.py.ipynb"

rule tsvs2hdf5:
	input:
		rules.vienna_cmbn.output,
		"results/filter/main.tsv",
		"results/extened/lineage.tsv"
	output:
		"results/final/result.hdf5"
	threads:
		1
	conda:
		"../envs/add_info.yaml"
	notebook:
		"../notebooks/tsv2hdf5.py.ipynb"
