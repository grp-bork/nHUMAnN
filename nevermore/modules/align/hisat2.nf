process hisat2_build {
    container "quay.io/biocontainers/hisat2:2.2.1--hdbdd923_6"
    // we need a hisat2/samtools mixed container
    // container "registry.git.embl.de/schudoma/hisat2-docker:latest"

    input:
    tuple val(sample), path(genomeseq)

    output:
    tuple val(sample), path("${sample.id}/hisat2/${sample.id}*"), emit: index

    script:
    """
    mkdir -p ${sample.id}/hisat2/

    gzip -dc ${genomeseq} > genome.fa

    hisat2-build -f genome.fa ${sample.id}/hisat2/${sample.id}

    rm -vf genome.fa
    """

}

process hisat2_align {
    // container "quay.io/biocontainers/hisat2:2.2.1--hdbdd923_6"
    // we need a hisat2/samtools mixed container
    container "registry.git.embl.de/schudoma/hisat2-docker:latest"

    input:
    tuple val(sample), path(fastqs), path(index)

    output:
    // tuple val(sample), path("${sample.id}/hisat2_align/${sample.id}.bam"), path("${sample.id}/hisat2_align/${sample.id}.bam.bai"), emit: bam
    tuple val(sample), path("${sample.id}/hisat2_align/${sample.id}.bam"), emit: bam
    tuple val(sample), path("${sample.id}/hisat2_align/${sample.id}.bam.bai"), emit: bai
    tuple val(sample), path("${sample.id}.HISAT2.DONE"), emit: sentinel

    script:

    def input_files = ""
	def r1_files = fastqs.findAll( { it.name.endsWith("_R1.fastq.gz") && !it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )
	def r2_files = fastqs.findAll( { it.name.endsWith("_R2.fastq.gz") } )
	def orphans = fastqs.findAll( { it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )

	if (r1_files.size() != 0 && r2_files.size() != 0) {
		input_files += "-1 ${r1_files.join(' ')} -2 ${r2_files.join(' ')}"
		single_reads = false
	} else if (r1_files.size() != 0) {
		input_files += "-U ${r1_files.join(' ')}"
	} else if (r2_files.size() != 0) {
		input_files += "-U ${r2_files.join(' ')}"
	} else if (orphans.size() != 0) {
		input_files += "-U ${orphans.join(' ')}"
	}

    // --fr/--rf/--ff
    def threads = task.cpus.intdiv(2)
    def hisat2_options = "-p ${threads} -q --phred33"
    if (params.hisat2_no_spliced_alignment) {
        // --no-spliced-alignment
        hisat2_options += " --no-spliced-alignment"
    }

    def index_id = index[0].name.replaceAll(/.[0-9]+.ht2$/, "")
    // -S ${sample.id}/hisat2_align/${sample.id}.sam
    // index_id=\$(ls ${index[0]} | sed 's/\\.[0-9]\\+\\.ht2\$//')
    """
    mkdir -p ${sample.id}/hisat2_align/ tmp/

    export TMPDIR=tmp/

    hisat2 -x ${index_id} ${hisat2_options} ${input_files} > ${sample.id}.sam
    samtools sort -@ ${threads} ${sample.id}.sam > ${sample.id}/hisat2_align/${sample.id}.bam
    samtools index ${sample.id}/hisat2_align/${sample.id}.bam
    rm -fv ${sample.id}.sam

    touch ${sample.id}.HISAT2.DONE
    """
    // echo "tmpdir is \$TMPDIR"

    // hisat2 -x ${sample.id} ${hisat2_options} ${input_files} > ${sample.id}.sam
}



