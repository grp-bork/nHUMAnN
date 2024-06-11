process fq2fa {
	container "quay.io/biocontainers/bbmap:39.06--h92535d8_0"

    input:
    tuple val(sample), path(fq)

    output:
    tuple val(sample), path("out/${sample.id}*.fasta"), emit: reads

    script:
	def maxmem = task.memory.toGiga()
	def r2 = (sample.is_paired) ? "in2=${sample.id}_R2.fastq.gz out2=out/${sample.id}_R2.fasta" : ""
	def qual_modifier = ""
	if (params.pb_reads) {
		qual_modifier = "qin=33"
	}

	"""
	set -e -o pipefail
	mkdir -p out/
	reformat.sh -Xmx${maxmem}g in=${sample.id}_R1.fastq.gz out=out/${sample.id}_R1.fasta ${r2} trimreaddescription=t ${qual_modifier}
	"""
	// https://forum.qiime2.org/t/bug-q2-itsxpresss-dependency-bbmap-cannot-handle-pacbio-ccs/17612/4
}
