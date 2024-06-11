process fq2bam {
	container "quay.io/biocontainers/bbmap:39.06--h92535d8_0"

    input:
    tuple val(sample), path(fq)

    output:
    tuple val(sample), path("out/${sample.id}.bam"), emit: reads

    script:
	def maxmem = task.memory.toGiga()
	def r2 = (sample.is_paired) ? "in2=${sample.id}_R2.fastq.gz" : ""
	def qual_modifier = ""
	if (params.pb_reads) {
		qual_modifier = "qin=33"
	}

	"""
	set -e -o pipefail
	mkdir -p out/
	reformat.sh -Xmx${maxmem}g in=${sample.id}_R1.fastq.gz ${r2} trimreaddescription=t out=stdout.bam ${qual_modifier} | samtools addreplacerg -r "ID:A" -r "SM:${sample.id}" -o out/${sample.id}.bam -
	"""
	// https://forum.qiime2.org/t/bug-q2-itsxpresss-dependency-bbmap-cannot-handle-pacbio-ccs/17612/4
}
