process prepare_fastqs {
    input:
    tuple val(sample), path(fq)

    output:
    tuple val(sample), path("fastq/${sample.id}/${sample.id}_R*.fastq.gz"), emit: reads

    script:
    if (sample.is_paired) {
        """
        mkdir -p fastq/${sample.id}
        ln -sf ../../${fq[0]} fastq/${sample.id}/${sample.id}_R1.fastq.gz
        ln -sf ../../${fq[1]} fastq/${sample.id}/${sample.id}_R2.fastq.gz
        """
    } else {
        """
        mkdir -p fastq/${sample.id}
        ln -sf ../../${fq[0]} fastq/${sample.id}/${sample.id}_R1.fastq.gz
        """
    }
}
