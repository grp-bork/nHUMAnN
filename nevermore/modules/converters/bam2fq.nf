process bam2fq {
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"
    input:
    tuple val(sample), path(bam)
    val(keep_unmapped)
    label "process_high"

    output:
    tuple val(sample), path("fastq/${sample.id}/${sample.id}*.fastq.gz"), emit: reads

    script:

    def filter_flags = "-F 0x900"
    if (keep_unmapped == true) {
        if (sample.is_paired) {
            filter_flags = "-F 0x900 -f 0xc"
        } else {
            filter_flags = "-F 0x900 -f 0x4"
        }
    }

    """
    set -o pipefail
    mkdir -p fastq/${sample.id} tmp/
    samtools collate -T tmp/tmpfile -@ $task.cpus -u -O $bam | samtools fastq ${filter_flags} -0 ${sample.id}_other.fastq.gz -1 ${sample.id}_R1.fastq.gz -2 ${sample.id}_R2.fastq.gz

    if [[ "\$?" -eq 0 ]];
    then

        if [[ -z "\$(gzip -dc ${sample.id}_R1.fastq.gz | head -n 1)" ]];
        then
            if [[ ! -z "\$(gzip -dc ${sample.id}_other.fastq.gz | head -n 1)" ]];
            then
                mv -v ${sample.id}_other.fastq.gz fastq/${sample.id}/${sample.id}_R1.fastq.gz;
            fi;
        else
                mv -v ${sample.id}_R1.fastq.gz fastq/${sample.id}/;
                if [[ ! -z "\$(gzip -dc ${sample.id}_R2.fastq.gz | head -n 1)" ]];
                then
                    mv -v ${sample.id}_R2.fastq.gz fastq/${sample.id}/;
                fi;
        fi;

        ls -l *.fastq.gz
        ls -l fastq/${sample.id}/*.fastq.gz
        rm -rf *.fastq.gz
    fi;
    """
}
