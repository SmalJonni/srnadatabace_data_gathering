#   ___ ___ _____   ___   _ _____ _   
#  / __| __|_   _| |   \ /_\_   _/_\  
# | (_ | _|  | |   | |) / _ \| |/ _ \ 
#  \___|___| |_|   |___/_/ \_\_/_/ \_\
import csv
localrules: get_SRA_sample_list, retrieve_data, get_reference_genomes

#GET SRA DATA

rule get_SRA_sample_list:
	output:
		"results/samples/list.txt"
	conda:
		"../envs/denovo_preprocess.yaml"
	threads:
		1
	params:
		address=config["email_address"]
	notebook:
		"../notebooks/get_data_list.py.ipynb"

#quite ineffective since this function is called multiple times
sample_lookup={}
def get_table():
	in_file=checkpoints.get_reference_genomes.get().output[0]
	in_file=rules.get_SRA_sample_list.output[0]
	global sample_lookup
	if not sample_lookup:
		with open(in_file, newline='') as csvfile:
			spamreader = csv.DictReader(csvfile, delimiter='\t')
			for row in spamreader:
				sample_lookup[row["Accession"]]=row["TaxID"]
	return sample_lookup

#GET DATA
checkpoint rule retrieve_data:
	input:
		rules.get_SRA_sample_list.output
	output:
		directory("results/data/"),
		temp(directory("results/prefetches/"))
	conda:
		"../envs/denovo_preprocess.yaml"
	threads:
		1
	shadow:
		"shallow"
	params:
		address=config["email_address"]
	shell:
		"""
		mkdir -p {output}
		tail -n+2 {input} | awk '{{print $1}}' > tmp_prefetch
		prefetch --option-file tmp_prefetch   --output-directory {output[1]}
		fasterq-dump --threads {threads}  --skip-technical   --split-3   -O {output[0]}   {output[1]}/*
		"""

#GET REF GENOMES & gbffs
checkpoint rule get_reference_genomes:
	input:
		#This dependency is artificially added to make sure that not too many threads are querying from NCBI at the same time and lead to timeouts 
		rules.get_SRA_sample_list.output,
		rules.retrieve_data.output,
		rules.prepare_pytaxonkit_place.output, #Note this comes from another file
	output:
		directory("results/ref_genomes/"),
		directory("results/gbffs/")
	conda:
		"../envs/denovo_preprocess.yaml"
	threads:
		1
	params:
		address=config["email_address"]
	notebook:
		"../notebooks/get_ref_genomes.py.ipynb"

#TODO come up with a more elegenant solution		
rule get_longest_entry:
	input:
		"results/ref_genomes/{wc}/genome.fa.gz"
	output:
		"results/ref_genomes_longest/{wc}/genome.fa"
	conda:
		"../envs/denovo_preprocess.yaml"
	threads:
		1
	notebook:
		"../notebooks/extract_long.py.ipynb"


#TODO come up with a more elegenant solution		
rule get_reference_genome_size:
	input:
		"results/ref_genomes_longest/{wc}/genome.fa"
	output:
		"results/ref_genomes_longest/stat/{wc}.txt"
	conda:
		"../envs/denovo_preprocess.yaml"
	threads:
		1
	shell:
		"cat {input}|grep -v '>' | tr -d '[:space:]' | wc -c > {output} "
		
#   ___   ___ 
#  / _ \ / __|
# | (_) | (__ 
#  \__\_\\___|

#QC READS
#Remove adapter with trimgalore
#TODO use inheritance instead of copy & pasting for pe&se ..
rule trim_galore_pe:
	input:
		["results/data/{sample}_1.fastq", "results/data/{sample}_2.fastq"]
	output:
		"results/trimmed/{sample}_1_val_1.fq",
		"results/trimmed/{sample}_2_val_2.fq",
		"results/trimmed/{sample}_1.fastq_trimming_report.txt",
		"results/trimmed/{sample}_2.fastq_trimming_report.txt",
	params:
		extra="--small_rna",
	log:
		"logs/trim_galore/{sample}.log"
	threads:
		8
	wrapper:
		"v1.7.0/bio/trim_galore/pe"

rule trim_galore_se:
	input:
		"results/data/{sample}.fastq",
	output:
		"results/trimmed/{sample}_trimmed.fq",
		"results/trimmed/{sample}.fastq_trimming_report.txt",
	params:
		extra="--small_rna",
	log:
		"logs/trim_galore/{sample}.log"
	threads:
		8
	wrapper:
		"v1.7.0/bio/trim_galore/se"


#    _   _ _                         _   
#   /_\ | (_)__ _ _ _  _ __  ___ _ _| |_ 
#  / _ \| | / _` | ' \| '  \/ -_) ' \  _|
# /_/ \_\_|_\__, |_||_|_|_|_\___|_||_\__|
#           |___/                        

rule bowtie2_build:
	input:
		ref="results/ref_genomes_longest/{wc}/genome.fa",
	output:
		multiext(
			"results/ref_genomes_longest/{wc}/genome.fa",
			".1.bt2",
			".2.bt2",
			".3.bt2",
			".4.bt2",
			".rev.1.bt2",
			".rev.2.bt2",
		),
	log:
		"logs/bowtie2_build/build_{wc}.log",
	params:
		extra="",  # optional parameters
	threads: 
		10
	wrapper:
		"v1.7.0/bio/bowtie2/build"

rule bowtie2_pe:
	input:
		sample=rules.trim_galore_pe.output[0:1],
		idx=lambda wildcards: multiext("results/ref_genomes_longest/"+str(get_table()[wildcards.sample])+"/genome.fa",
			".1.bt2",
			".2.bt2",
			".3.bt2",
			".4.bt2",
			".rev.1.bt2",
			".rev.2.bt2",
			)
	output:
		"results/alignment/{sample}.bam"
	log:
		"logs/alignment/{sample}.bam"
	params:
		extra=""  # optional parameters
	threads: 20  # Use at least two threads
	wrapper:
		"v1.7.0/bio/bowtie2/align"

use rule bowtie2_pe as bowtie2_se with:
	input:
		sample=[rules.trim_galore_se.output[0]],
		idx=lambda wildcards: multiext("results/ref_genomes_longest/"+str(get_table()[wildcards.sample])+"/genome.fa",
			".1.bt2",
			".2.bt2",
			".3.bt2",
			".4.bt2",
			".rev.1.bt2",
			".rev.2.bt2",
			)

#  ___            _ _    _   _          
# | _ \_ _ ___ __| (_)__| |_(_)___ _ _  
# |  _/ '_/ -_) _` | / _|  _| / _ \ ' \ 
# |_| |_| \___\__,_|_\__|\__|_\___/_||_|
			                           
#Convert gbff

rule xtractgbff:
	input:
		gbff="results/gbffs/{wc}.gbff.gz",
		fa=rules.get_longest_entry.output
	output:
		"results/gbffs_long/{wc}.gbff"
	threads:
		1
	shell:
		"""
		mkdir -p $(dirname {output})
		LINE=$(zcat  {input.gbff} | grep -n LOCUS | sort -n -k3 | tail -n1 | awk -F":" '{{print $1}}')
		ELINE=$(zcat {input.gbff}  | tail -n+$LINE  | grep -n //$ | awk -F":" '{{print $1}}' | head -n1)
		head -n $ELINE <(zcat {input.gbff}  | tail -n+$LINE)  >  {output}
		"""
#TODO improve this terrible script above .. .
rule gbff2ppt:
	input:
		rules.xtractgbff.output
	output:
		"results/ptts/{wc}.ptt"
	conda:
		"../envs/gbff2ptt.yml"
	threads:
		1
	shell:
		"""
		perl workflow/scripts/gb2ptt.pl --infile {input}
		mv {input}.ptt {output}
		"""
#TODO  do the above smarter


#INSTALLS APERO and creates a flag file
rule install_APERO:
	output:
		"results/flag/apero_installed.txt"
	conda:
		"../envs/denovo_preprocess.yaml"
	threads:
		1
	notebook:
		"../notebooks/install_apero.r.ipynb"

rule predict_novel_pe:
	input:
		ptt=lambda wildcards: "results/ptts/"+get_table()[wildcards.sample]+".ptt",
		bam=rules.bowtie2_pe.output,
		flag=rules.install_APERO.output,
		reads=rules.trim_galore_pe.output,
		genome_size=lambda wildcards: "results/ref_genomes_longest/stat/"+get_table()[wildcards.sample]+".txt"
	output:
		"results/APERO/{sample}.gff"
	conda:
		"../envs/denovo_preprocess.yaml"
	params:
		pe=True
	threads:
		20
#	shadow:
#		"shallow"
	notebook:
		"../notebooks/predict_apero.r.ipynb"

use rule predict_novel_pe as predict_novel_se with:
	input:
		ptt=lambda wildcards: "results/ptts/"+get_table()[wildcards.sample]+".ptt",
		bam=rules.bowtie2_se.output,
		flag=rules.install_APERO.output,
		reads=rules.trim_galore_se.output,
		genome_size=lambda wildcards: "results/ref_genomes_longest/stat/"+get_table()[wildcards.sample]+".txt"
	params:
		pe=False

rule count_reads_pe:
	input:
		reads=rules.trim_galore_pe.output
	output:
		read_count="results/read_count/{sample}.txt"
	shell:
		"""
		A=$( cat  {input[0]} | wc -l  )
		echo $((A/2)) >{output}
		"""

rule count_reads_se:
	input:
		reads=rules.trim_galore_se.output
	output:
		read_count="results/read_count/{sample}.txt"
	shell:
		"""
		A=$( cat  {input[0]} | wc -l  )
		echo $((A/4)) >{output}
		"""

rule filter_novel_predictions:
	input:
		apero="results/APERO/{sample}.gff",
		genome_size=lambda wildcards: "results/ref_genomes_longest/stat/"+get_table()[wildcards.sample]+".txt",
		read_count="results/read_count/{sample}.txt"
	output:
		"results/APERO_filtered/{sample}.gff"
	conda:
		"../envs/denovo_preprocess.yaml"
	threads:
		1
	notebook:
		"../notebooks/filter_predictions.r.ipynb"


rule extract_sequence:
	input:
		gff=rules.filter_novel_predictions.output,
		FASTA=lambda wildcards: "results/ref_genomes_longest/"+str(get_table()[wildcards.sample])+"/genome.fa"
	output:
		"results/APERO_sequences/{sample}.fa"
	conda:
		"../envs/denovo_preprocess.yaml"
	threads:
		1
	shadow:
		"shallow"
	shell:
		"""
		CHROM=$(head -n1 {input.FASTA} |  awk -F" " '{{print $1}}' | sed 's/^>//g')
		tail -n+2 {input.gff}    | awk -v a=$CHROM '{{print a"\\t"$2"\\t"($2+$4)}}' | sed 's/\\.5\\t/\\t/g' | sed 's/\\.5\$//g'   > tmp.bed
		bedtools getfasta -fi {input.FASTA} -bed tmp.bed > {output}
		"""
#No idea but apero thinks it is funny to include .5 positions??
rule prepare_denovo:
	input:
		rules.extract_sequence.output,
		gbgff=lambda wildcards: "results/gbffs_long/"+str(get_table()[wildcards.sample])+".gbff",
		pytaxonkit=rules.prepare_pytaxonkit_place.output #Note this comes from another file
	output:
		"results/denovo_prepared/{sample}.tsv"
	conda:
		"../envs/denovo_preprocess.yaml"
	threads:
		1
	notebook:
		"../notebooks/prepare_denovo.py.ipynb"


rule combine_predictions:
	input:
		lambda wildcard: expand("results/denovo_prepared/{sample}.tsv",sample=get_table().keys())
	output:
		"result/cmbn/denovo.tsv"
	conda:
		"../envs/denovo_preprocess.yaml"
	threads:
		1
	notebook:
		"../notebooks/cbmn_denovo.r.ipynb"

