process minimap2_align {
	container "registry.git.embl.de/schudoma/minimap2-docker:latest"
	label 'align'

	input:
	tuple val(sample), path(fastqs)
	path(reference)
	val(do_name_sort)

	output:
	tuple val(sample), path("${sample.id}/${sample.id}.sam"), emit: sam

	script:
	// def reads = (sample.is_paired) ? "${sample.id}_R1.fastq.gz ${sample.id}_R2.fastq.gz" : "${sample.id}_R1.fastq.gz"

	def input_files = ""
	def r1_files = fastqs.findAll( { it.name.endsWith("_R1.fastq.gz") && !it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )
	def r2_files = fastqs.findAll( { it.name.endsWith("_R2.fastq.gz") } )
	def orphans = fastqs.findAll( { it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )

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




	def threads = task.cpus.intdiv(2)
	// def mm_options = "--sam-hit-only -t ${threads} -x sr --secondary=yes -a"
	def mm_options = "-t ${threads} -x sr --secondary=yes -a"

	// def sort_cmd = "| " + ((do_name_sort) ? "samtools collate -@ ${threads} -o ${sample.id}.bam - tmp/collated_bam" : "samtools sort -@ ${threads} -o ${sample.id}.bam -")
	def sort_cmd = ""  // we cannot convert large catalogue alignments to bam, hence we cannot properly sort those

	"""
	set -e -o pipefail

	mkdir -p ${sample.id}/ tmp/
	minimap2 ${mm_options} --split-prefix ${sample.id}_split ${reference} ${input_files} ${sort_cmd} > ${sample.id}/${sample.id}.sam

	rm -rvf tmp/
	"""
}