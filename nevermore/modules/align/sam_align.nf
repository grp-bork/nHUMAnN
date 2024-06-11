process minimap2_align {
	container "quay.io/biocontainers/minimap2:2.28--he4a0461_0"
	label 'align'

	input:
	tuple val(sample), path(reads)
	path(reference)

	output:
	tuple val(sample), path("${sample.id}/${sample.id}.sam"), emit: sam

	script:
	def reads = (sample.is_paired) ? "${sample.id}_R1.fastq.gz ${sample.id}_R2.fastq.gz" : "${sample.id}_R1.fastq.gz"
	def mm_options = "--sam-hit-only -t ${task.cpus} -x sr --secondary=yes -a"

	"""
	mkdir -p ${sample.id}/
	minimap2 ${mm_options} --split-prefix ${sample.id}_split ${reference} ${reads} > ${sample.id}/${sample.id}.sam
	"""
}


process bwa_mem_align {
	container "quay.io/biocontainers/bwa:0.7.3a--he4a0461_9"
	label 'align'

	input:
	tuple val(sample), path(reads)
	path(reference)

	output:
	tuple val(sample), path("${sample.id}/${sample.id}.sam"), emit: sam

	script:
	def reads = (sample.is_paired) ? "${sample.id}_R1.fastq.gz ${sample.id}_R2.fastq.gz" : "${sample.id}_R1.fastq.gz"
	def bwa_options = "-a -t ${task.cpus} -K 10000000"

	"""
	mkdir -p ${sample.id}/
	bwa mem ${bwa_options} \$(readlink ${reference}) ${reads} |\
		awk -F '\t' -v OFS='\t' 'substr(\$1,1,1) == "@" || !and(\$2, 4)' > ${sample.id}/${sample.id}.sam
	"""

}
