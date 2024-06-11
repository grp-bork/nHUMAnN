process run_metaphlan3 {
	container "quay.io/biocontainers/metaphlan:3.1.0--pyhb7b1952_0"
	
	input:
	tuple val(sample), path(fastqs)
	path(mp3_db)

	output:
	tuple val(sample), path("metaphlan3_tables/${sample.id}.mp3.txt"), emit: mp3_table
	// tuple val(sample), path("${sample.id}.bowtie2.bz2"), emit: mp3_bt2
	
	script:
	def mp3_params = "--index ${params.mp3_db_version} --bowtie2db \$(readlink ${mp3_db}) --input_type fastq --nproc ${task.cpus} --tmp_dir tmp/"
	def mp3_input = ""
	def bt2_out = "--bowtie2out ${sample.id}.bowtie2.bz2"

	
	if (fastqs instanceof Collection && fastqs.size() >= 2) {
		mp3_input = "${sample.id}_R1.fastq.gz,${sample.id}_R2.fastq.gz"
	// } else if (fastqs instanceof Collection && fastqs.size() == 3) {
	// 	mp3_input = "${sample.id}_R1.fastq.gz,${sample.id}_R2.fastq.gz,${sample.id}.singles_R1.fastq.gz"
	} else {
		mp3_input = "${fastqs}"
	}

	def additional_mp3_params = ""
	if (params.mp3_params) {
		additional_mp3_params = params.mp3_params
	}


	"""
	mkdir -p tmp/ metaphlan3_tables/

	final_input=${mp3_input}

	if [[ -e ${sample.id}.singles_R1.fastq.gz && ${mp3_input} != ${sample.id}.singles_R1.fastq.gz ]]; then
		nlines=\$(gzip -dc ${sample.id}.singles_R1.fastq.gz | head | wc -l)
		if [[ \$nlines -gt 4 ]]; then
			final_input=\$final_input\",${sample.id}.singles_R1.fastq.gz\"
		fi
	fi

	echo metaphlan \$final_input ${additional_mp3_params} ${mp3_params} ${bt2_out} -o metaphlan3_tables/${sample.id}.mp3.txt

	metaphlan \$final_input ${additional_mp3_params} ${mp3_params} ${bt2_out} -o metaphlan3_tables/${sample.id}.mp3.txt
	"""
}

process combine_metaphlan3 {
	container "quay.io/biocontainers/metaphlan:3.1.0--pyhb7b1952_0"

	input:
	tuple val(sample), path(bt2)

	output:
	tuple val(sample), path("metaphlan3/${sample.id}/${sample.id}.mp3.txt"), emit: mp3_table

	script:
	def mp3_params = "--input_type bowtie2out --nproc ${task.cpus} --tmp_dir tmp/"
	def bt2_out = "--bowtie2out ${sample.id}.bowtie2.bz2"
	def mp3_input = "${sample.id}.bowtie2.bz2,${sample.id}.singles.bowtie2.bz2"
	"""
	mkdir -p metaphlan3/${sample.id}/

	metaphlan ${mp3_input} ${mp3_params} -o metaphlan3/${sample.id}/${sample.id}.mp3.txt
	"""
}


process collate_metaphlan3_tables {
	container "quay.io/biocontainers/metaphlan:3.1.0--pyhb7b1952_0"

	input:
	path(tables)

	output:
	path("metaphlan3_abundance_table.txt")

	script:
	"""
	mkdir -p metaphlan3/

	merge_metaphlan_tables.py ${tables} > metaphlan3_abundance_table.txt
	"""

}