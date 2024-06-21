process kallisto_index {
	container "quay.io/biocontainers/kallisto:0.50.1--hc877fd6_1"
	label "medium"

	input:
	tuple val(sample), path(genes)

	output:
	tuple val(sample), path("kallisto/index/${sample.id}.idx"), emit: index

	script:
	"""
	mkdir -p kallisto/index/${sample.id}/

	kallisto index -i kallisto/index/${sample.id}.idx ${genes}
	"""
	
}

params.profilers = [:]
params.profilers.kallisto = [:]
params.profilers.kallisto.bootstrap = 100

process kallisto_quant {
	container "quay.io/biocontainers/kallisto:0.50.1--hc877fd6_1"
	label "medium"

	input:
	tuple val(sample), path(fastqs), path(kallisto_index)

	output:
	tuple val(sample), path("kallisto/quant/${sample.id}/*"), emit: quant

	script:

	def input_files = ""
	def r1_files = fastqs.findAll( { it.name.endsWith("_R1.fastq.gz") && !it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )
	def r2_files = fastqs.findAll( { it.name.endsWith("_R2.fastq.gz") } )
	def orphans = fastqs.findAll( { it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )

	def single_reads = true
	if (r1_files.size() != 0 && r2_files.size() != 0) {
		input_files += "${r1_files.join(' ')} ${r2_files.join(' ')}"
		single_reads = false
	} else if (r1_files.size() != 0) {
		input_files += "${r1_files.join(' ')}"
	} else if (r2_files.size() != 0) {
		input_files += "${r2_files.join(' ')}"
	} else if (orphans.size() != 0) {
		input_files += "${orphans.join(' ')}"
	}

	def single_flags = ""
	def calc_libstats = ""
	if (single_reads) {

		calc_libstats += "mean_length=\$(gzip -dc ${input_files} | awk 'NR%4==2 {sum_len+=length(\$0); n+=1; } END { printf('%d\n', sum_len / n); }')\n"
		calc_libstats += "std_length=\$(gzip -dc ${input_files} | awk -v mean=\$mean_length 'NR%4==2 {sum_len+=(length(\$0)-mean)**2; n+=1; } END { printf('%d\n', sqrt(sum_len/(n-1)))}')"
		single_flags += "--single"
		single_flags += " -l \$mean_length -s \$std_length"

	}

	"""
	mkdir -p kallisto/quant/${sample.id}/

	${calc_libstats}

	kallisto quant -i ${kallisto_index} -b ${params.profilers.kallisto.bootstrap} -o kallisto/quant/${sample.id} ${single_flags} ${input_files}
	"""	
	// salmon quant -p ${task.cpus} -i ${salmon_index} -l ${params.profilers.salmon.quant.libtype} ${input_files} --validateMappings -o salmon/quant/${sample.id}/
	// ./bin/salmon quant -i transcripts_index -l <LIBTYPE> -1 reads1.fq -2 reads2.fq --validateMappings -o transcripts_quant

}
