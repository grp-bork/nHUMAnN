process fq2bam {
	container "quay.io/biocontainers/bbmap:39.06--h92535d8_0"
	label "process_high"

    input:
    tuple val(sample), path(fastqs)

    output:
    tuple val(sample), path("out/${sample.id}.bam"), emit: reads

    script:
	def maxmem = task.memory.toGiga()
	// def r2 = (sample.is_paired) ? "in2=${sample.id}_R2.fastq.gz" : ""
	def qual_modifier = ""
	if (params.pb_reads) {
		qual_modifier = "qin=33"
	}

	def r1_files = fastqs.findAll( { it.name.endsWith("_R1.fastq.gz") && !it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )
	def r2_files = fastqs.findAll( { it.name.endsWith("_R2.fastq.gz") } )
	def orphans = fastqs.findAll( { it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )

	def input_files = ""
	if (r1_files.size() != 0 && r2_files.size() != 0) {
		// input_files += "-f ${r1_files.join(' ')} -r ${r2_files.join(' ')}"
		input_files = "in=${r1_files.join(',')} in2=${r2_files.join(',')}"
	} else if (r1_files.size() != 0) {
		// input_files += "-s ${r1_files.join(' ')}"
		input_files += "in=${r1_files.join(',')}"
	} else if (r2_files.size() != 0) {
		// input_files += "-s ${r2_files.join(' ')}"
		input_files += "in=${r2_files.join(',')}"
	} else if (orphans.size() != 0) {
		input_files += "in=${orphans.join(',')}"
	}


	"""
	set -e -o pipefail
	mkdir -p out/
	reformat.sh -Xmx${maxmem}g ${input_files} trimreaddescription=t out=stdout.bam ${qual_modifier} | samtools addreplacerg -r "ID:A" -r "SM:${sample.id}" -o out/${sample.id}.bam -
	"""
	// reformat.sh -Xmx${maxmem}g in=${sample.id}_R1.fastq.gz ${r2} trimreaddescription=t out=stdout.bam ${qual_modifier} | samtools addreplacerg -r "ID:A" -r "SM:${sample.id}" -o out/${sample.id}.bam -
	// https://forum.qiime2.org/t/bug-q2-itsxpresss-dependency-bbmap-cannot-handle-pacbio-ccs/17612/4
}
