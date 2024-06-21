process bowtie2_build {
	container "quay.io/biocontainers/bowtie2:2.5.3--py39h6fed5c7_1"
	tag "${sample.id}"

	input:
    tuple val(sample), path(genomeseq)

    output:
    tuple val(sample), path("${sample.id}/bowtie2/${sample.id}*"), emit: index

    script:
    """
    mkdir -p ${sample.id}/bowtie2/

    gzip -dc ${genomeseq} > genome.fa

    bowtie2-build --threads ${task.cpus} -f genome.fa ${sample.id}/bowtie2/${sample.id}

    rm -vf genome.fa
    """

}


process bowtie2_align {
	// container "quay.io/biocontainers/bowtie2:2.5.3--py39h6fed5c7_1"
	container "registry.git.embl.de/schudoma/bowtie2-docker:latest"
	tag "${sample.id}"

	input:
    tuple val(sample), path(fastqs), path(index)

    output:
    tuple val(sample), path("${sample.id}/bowtie2_align/${sample.id}.bam"), emit: bam
    tuple val(sample), path("${sample.id}/bowtie2_align/${sample.id}.bam.bai"), emit: bai
    tuple val(sample), path("${sample.id}.BOWTIE2.DONE"), emit: sentinel

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
    def bowtie2_options = "-p ${threads} -q --phred33"
    
    def index_id = index[0].name.replaceAll(/.[0-9]+.bt2[l]?$/, "")
    // -S ${sample.id}/hisat2_align/${sample.id}.sam
    // index_id=\$(ls ${index[0]} | sed 's/\\.[0-9]\\+\\.ht2\$//')
    """
    mkdir -p ${sample.id}/bowtie2_align/ tmp/

    export TMPDIR=tmp/

    bowtie2 -x ${index_id} ${bowtie2_options} ${input_files} > ${sample.id}.sam
    samtools sort -@ ${threads} ${sample.id}.sam > ${sample.id}/bowtie2_align/${sample.id}.bam
    samtools index ${sample.id}/bowtie2_align/${sample.id}.bam
    rm -fv ${sample.id}.sam

    touch ${sample.id}.BOWTIE2.DONE
    """

}