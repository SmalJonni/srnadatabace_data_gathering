
rule get_sRNATarbase_info:
	output:
		"results/Database_knowledge/sRNAtarbase.tsv"
	conda:
		"../envs/query_databases.yaml"
	log:
		notebook="logs/notebooks/get_sRNAtarbase.py.ipynb"
	threads:
		1
	notebook:
		"../notebooks/get_sRNAtarbase_processed.py.ipynb"

rule get_RNA_central:
	output:
		"results/Database_knowledge/RNAcentral_raw.tsv"
	conda:
		"../envs/query_databases.yaml"
	threads:
		1
	shell:
		"""
		export PGPASSWORD=NWDMCE5xdipIjRrp
		psql -h hh-pgsql-public.ebi.ac.uk  -p 5432 -d pfmegrnargs -U reader    -f  resources/rnacentral_query.sql -o {output}
		"""
rule mod_RNA_central:
	input:
		"results/Database_knowledge/RNAcentral_raw.tsv"
	output:
		"results/Database_knowledge/RNAcentral.tsv"
	conda:
		"../envs/query_databases.yaml"

	threads:
		1
	shell:
		"""
		echo "OSID\tTaxa\tSequence\tFunction\tSource" > {output}
		tail -n+3 {input} | head -n-1 | sed 's/|/\t/g' | sed 's/ //g'| awk '{{print $0"\tRNAcentral"}}' >> {output}
		"""


rule get_Gene_DB:
	output:
		"results/Database_knowledge/GeneDB.tsv"
	params:
		address=config["email_address"]
	conda:
		"../envs/query_databases.yaml"
	log:
		notebook="logs/notebooks/get_gene_db.py_processed.ipynb"
	threads:
		1
	notebook:
		"../notebooks/get_gene_db.py.ipynb"

rule get_CoryneRegNet:
	output:
		"results/Database_knowledge/CoryneRegNet.tsv"
	shadow:
		"minimal"
	threads:
		1
	shell:
		"""
		wget -O cory.zip  https://exbio.wzw.tum.de/coryneregnet/downloadFiles.htm?fileName=AllOrganismsFiles.zip
		unzip cory.zip
		echo  'Taxa\tSequence\tSource\tFunction' > {output}
		cat *_rna.csv | sed '/^locus_tag/d'   | awk -F"\t" '{{print $6"\t"$11"\tCorynet\tidk"}}' >> {output}  
		"""

rule get_Regulon_DB:
	output:
		"results/Database_knowledge/Regulondb.tsv"
	conda:
		"../envs/query_databases.yaml"
	log:
		notebook="logs/notebooks/get_regulondb_processed.py.ipynb"
	threads:
		1
	notebook:
		"../notebooks/get_regulondb.py.ipynb"


rule combine_database_knowledge:
	input:
		rules.get_sRNATarbase_info.output,
		rules.mod_RNA_central.output,
		#rules.get_Gene_DB.output,
		rules.get_CoryneRegNet.output,
		rules.get_Regulon_DB.output
	output:
		"results/Database_knowledge/combined.tsv"
	conda:
		"../envs/query_databases.yaml"
	log:
		notebook="logs/notebooks/combine_databases_processed.py.ipynb"
	threads:
		1
	notebook:
		"../notebooks/combine_databases.py.ipynb"

rule get_taxonomy_ontology:
	output:
		"results/Database_knowledge/ncbitaxon.owl"
	threads:
		1
	shell:
		"wget -O {output} http://purl.obolibrary.org/obo/ncbitaxon.owl"
