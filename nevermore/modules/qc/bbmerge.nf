process qc_bbmerge {
    container "quay.io/biocontainers/bbmap:39.06--h92535d8_0"
    label "bbduk"

    input:
    tuple val(sample), path(reads)

    output:
    tuple val("${sample}.singles"), path("${sample}.singles/${sample}.singles_M.fastq.gz"), optional: true, emit: merged
    tuple val(sample), path("${sample}/${sample}_R*.fastq.gz"), optional: true, emit: pairs

    script:
    def maxmem = task.memory.toGiga()
    def merge_params = "rsem=t extend2=20 iterations=5 ecct vstrict"

    """
    mkdir -p ${sample} ${sample}.merged

    bbmerge.sh -Xmx${maxmem}g t=${task.cpus} ${merge_params} in=${sample}_R1.fastq.gz in2=${sample}_R2.fastq.gz out=${sample}.singles/${sample}.singles_M.fastq.gz outu1=${sample}/${sample}_R1.fastq.gz outu2=${sample}/${sample}_R2.fastq.gz
    """
}

params.bbmerge_insert_size = "ecct extend2=20 iterations=5 k=62 adapter=default"


process qc_bbmerge_insert_size {
    container "quay.io/biocontainers/bbmap:39.06--h92535d8_0"
    label "bbduk"

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample.id}/bbmerge/${sample.id}.ihist.txt"), emit: isize_hist
    tuple val(sample), path("${sample.id}/bbmerge/${sample.id}.inserts.txt"), emit: inserts
    tuple val(sample), path("${sample.id}/bbmerge/${sample.id}.adapters.txt"), emit: adapters

    script:
    def maxmem = task.memory.toGiga()
    def compression = (reads[0].name.endsWith("gz")) ? "gz" : "bz2"

    def r1_files = reads.findAll( { it.name.endsWith("_R1.fastq.${compression}") } )
	def r2_files = reads.findAll( { it.name.endsWith("_R2.fastq.${compression}") } )

    def read1 = ""
    def read2 = ""
    def orphans = ""
    if (r1_files.size() != 0) {
        read1 += "in1=${r1_files[0]}"
        if (r2_files.size() != 0) {
            read2 += "in2=${r2_files[0]}"        
        }
    }

    """
    mkdir -p ${sample.id}/bbmerge/

    bbmerge.sh -Xmx${maxmem}g t=${task.cpus} ${read1} ${read2} \
      ihist=${sample.id}/bbmerge/${sample.id}.ihist.txt \
      outinsert=${sample.id}/bbmerge/${sample.id}.inserts.txt \
      outadapter=${sample.id}/bbmerge/${sample.id}.adapters.txt \
      ${params.bbmerge_insert_size}
    """
}