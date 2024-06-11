process merge_single_fastqs {
    container "quay.io/biocontainers/bbmap:39.06--h92535d8_0"

    input:
    tuple val(sample), path(fastqs)

    output:
    tuple val(sample), path("merged/${sample.id}_R1.fastq.gz"), emit: fastq

    script:

    def fastq_in = ""
    def prefix = ""
    def suffix = ""
    if (fastqs instanceof Collection && fastqs.size() == 2) {
        prefix = "cat ${fastqs} | gzip -dc - |"
        fastq_in = "stdin.fastq"
        suffix = " || { ec=\$?; [ \$ec -eq 141 ] && true || (exit \$ec); }"
    } else {
        fastq_in = "${fastqs[0]}"
    }

    """
    set -e -o pipefail
    mkdir -p merged/

    ${prefix} sortbyname.sh in=${fastq_in} out=merged/${sample.id}_R1.fastq.gz ${suffix}
    """
    // https://stackoverflow.com/questions/22464786/ignoring-bash-pipefail-for-error-code-141/72985727#72985727
}
