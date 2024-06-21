process flagstats {
    container "registry.git.embl.de/schudoma/align-docker:latest"
    label "default"

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample.id}/${sample.id}.flagstats.txt"), emit: flagstats
    tuple val(sample), path("${sample.id}/${sample.id}.libsize.txt"), emit: counts
    tuple val(sample), path("${sample.id}/${sample.id}.is_paired.txt"), emit: is_paired
    
    script:
    """
    mkdir -p ${sample.id}
    samtools flagstat $bam > "${sample.id}/${sample.id}.flagstats.txt"
    head -n 1 "${sample.id}/${sample.id}.flagstats.txt" | awk '{print \$1 + \$3}' > "${sample.id}/${sample.id}.libsize.txt"
    grep -m 1 "paired in sequencing" "${sample.id}/${sample.id}.flagstats.txt" | awk '{npaired = \$1 + \$3; if (npaired==0) {print "unpaired"} else {print "paired"};}' > "${sample.id}/${sample.id}.is_paired.txt"
    """
}
