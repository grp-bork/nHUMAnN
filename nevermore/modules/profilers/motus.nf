params.motus_tax_level = "mOTU"
params.motus_min_length = 75
params.motus_n_marker_genes = 3


process motus {
    container "quay.io/biocontainers/motus:3.1.0--pyhdfd78af_0"

    input:
    tuple val(sample), path(reads)
	path(motus_db)

    output:
    tuple val(sample), path("${sample.id}/${sample.id}.motus.txt"), emit: motus_out

    script:
    def motus_input = (sample.is_paired) ? "-f ${sample.id}_R1.fastq.gz -r ${sample.id}_R2.fastq.gz" : "-s ${sample.id}_R1.fastq.gz";
    """
    mkdir -p ${sample.id}
    motus profile -t $task.cpus -k ${params.motus_tax_level} -c -v 7 -q -l ${params.motus_min_length} -g ${params.motus_n_marker_genes} -db ${motus_db} ${motus_input} > ${sample.id}/${sample.id}.motus.txt
    """
}
