process qc_bbduk_stepwise_amplicon {
	container "quay.io/biocontainers/bbmap:39.06--h92535d8_0"
	label 'bbduk'

    input:
    tuple val(sample), path(reads)
	path(adapters)

    output:
    tuple val(sample), path("qc_reads/${sample.id}/${sample.id}_R*.fastq.gz"), emit: reads
    tuple val(sample), path("qc_reads/${sample.id}/${sample.id}.orphans_R1.fastq.gz"), emit: orphans, optional: true
    path("stats/qc/bbduk/${sample.id}/${sample.id}.*bbduk_stats.txt"), optional: true
    path("stats/qc/bbduk/${sample.id}/${sample.id}*lhist.txt"), emit: read_lengths, optional: true

    script:
	def maxmem = task.memory.toGiga()
	def compression = (reads[0].name.endsWith(".gz")) ? "gz" : "bz2"

	if (params.primers) {
		trim_params = "literal=${params.primers} minlength=${params.qc_minlen}"
	} else {
		trim_params = "ref=${adapters} minlength=${params.qc_minlen}"
	}

	def bbduk_call = "bbduk.sh -Xmx${maxmem}g t=${task.cpus} ordered=t trd=t"

	ref_p5_r1 = (params.primers) ? "literal=" + params.primers.split(",")[0] : "ref=${adapters}"
	ref_p5_r2 = (params.primers && !sample.is_paired) ? "literal=" + params.primers.split(",")[1] : "ref=${adapters}"
	ref_p3_r1 = ref_p5_r2
	ref_p3_r2 = ref_p5_r1

	def bbduk_full_call = ""
	def downstream_call = ""
	def lhist_call = "${bbduk_call} in1=qc_reads/${sample.id}/${sample.id}_R1.fastq.gz lhist=stats/qc/bbduk/${sample.id}/${sample.id}_R1.post_lhist.txt\n"

	if (sample.is_paired) {

		bbduk_full_call += "${bbduk_call} ${ref_p5_r1} minlength=${params.qc_minlen} ${params.p5_primer_params} in1=${sample.id}_R1.fastq.${compression} out1=fwd_p5.fastq.gz\n"
		bbduk_full_call += "${bbduk_call} ${ref_p5_r2} minlength=${params.qc_minlen} ${params.p5_primer_params} in1=${sample.id}_R2.fastq.${compression} out1=rev_p5.fastq.gz\n"

		if (params.long_reads) {

			bbduk_full_call += "${bbduk_call} ${ref_p3_r1} minlength=${params.qc_minlen} ${params.p3_primer_params} in1=fwd_p5.fastq.gz out1=fwd.fastq.gz\n"
			bbduk_full_call += "${bbduk_call} ${ref_p3_r2} minlength=${params.qc_minlen} ${params.p3_primer_params} in1=rev_p5.fastq.gz out1=rev.fastq.gz\n"

		} else {
			bbduk_full_call += "mv fwd_p5.fastq.gz fwd.fastq.gz\nmv rev_p5.fastq.gz rev.fastq.gz\n"
		}

		downstream_call += "gzip -dc fwd.fastq.gz | awk 'NR%4==1' | sed 's/^@//' | sed 's/\\/1//' | sort -T tmp/ > fwd.txt\n"
        downstream_call += "gzip -dc rev.fastq.gz | awk 'NR%4==1' | sed 's/^@//' | sed 's/\\/2//' | sort -T tmp/ > rev.txt\n"
		downstream_call += "join -1 1 -2 1 fwd.txt rev.txt > both.txt\n"
		downstream_call += "seqtk subseq fwd.fastq.gz both.txt | gzip -c - > qc_reads/${sample.id}/${sample.id}_R1.fastq.gz\n"
		downstream_call += "seqtk subseq rev.fastq.gz both.txt | gzip -c - > qc_reads/${sample.id}/${sample.id}_R2.fastq.gz\n"

		lhist_call += "${bbduk_call} in1=qc_reads/${sample.id}/${sample.id}_R2.fastq.gz lhist=stats/qc/bbduk/${sample.id}/${sample.id}_R2.post_lhist.txt\n"

	} else {
		bbduk_full_call += "${bbduk_call} ${trim_params} in1=${sample.id}_R1.fastq.${compression} out1=qc_reads/${sample.id}/${sample.id}_R1.fastq.gz stats=stats/qc/bbduk/${sample.id}/${sample.id}.fwd_bbduk_stats.txt lhist=stats/qc/bbduk/${sample.id}/${sample.id}.p5_lhist.txt\n"
	}


	"""
	mkdir -p ${sample.id}/
	mkdir -p stats/qc/bbduk/${sample.id}/
	mkdir -p qc_reads/${sample.id}/
	mkdir -p tmp/

	${bbduk_full_call}
	${downstream_call}
	${lhist_call}
	"""
}
