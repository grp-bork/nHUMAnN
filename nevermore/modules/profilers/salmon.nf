params.profilers.salmon.index.k = 31

process salmon_index {
	container "quay.io/biocontainers/salmon:1.10.3--hecfa306_0"

	input:
	tuple val(sample), path(genes)

	output:
	tuple val(sample), path("salmon/index/${sample.id}"), emit: index

	script:
	"""
	mkdir -p salmon/index/${sample.id}/

	salmon index -t ${genes} -i salmon/index/${sample.id}/ -k ${params.profilers.salmon.index.k}

	"""
	// > ./bin/salmon index -t transcripts.fa -i transcripts_index --decoys decoys.txt -k 31

}

params.profilers.salmon.quant.libtype = "IU"

process salmon_quant {
	container "quay.io/biocontainers/salmon:1.10.3--hecfa306_0"
	// label "align"

	input:
	tuple val(sample), path(fastqs), path(salmon_index)

	output:
	tuple val(sample), path("salmon/quant/${sample.id}/*"), emit: quant

	script:

	def input_files = ""
	def r1_files = fastqs.findAll( { it.name.endsWith("_R1.fastq.gz") && !it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )
	def r2_files = fastqs.findAll( { it.name.endsWith("_R2.fastq.gz") } )
	def orphans = fastqs.findAll( { it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )

	if (r1_files.size() != 0) {
		input_files += "-1 ${r1_files.join(' ')}"
	}
	if (r2_files.size() != 0) {
		input_files += " -2 ${r2_files.join(' ')}"
	}
	if (orphans.size() != 0) {
		input_files += " -r ${orphans.join(' ')}"
	}


	"""
	mkdir -p salmon/quant/${sample.id}/

	salmon quant -p ${task.cpus} -i ${salmon_index} -l ${params.profilers.salmon.quant.libtype} ${input_files} --validateMappings -o salmon/quant/${sample.id}/
	"""	
	// ./bin/salmon quant -i transcripts_index -l <LIBTYPE> -1 reads1.fq -2 reads2.fq --validateMappings -o transcripts_quant

}
