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
