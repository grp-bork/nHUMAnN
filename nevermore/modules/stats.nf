process flagstats {
    container "registry.git.embl.org/schudoma/align-docker:latest"
    label "default"

    input:
    tuple val(sample), path(bam)
    val(stage)

    output:
    tuple val(sample), path("${stage}/${sample.id}/${sample.id}.flagstats.txt"), emit: flagstats
    tuple val(sample), path("${stage}/${sample.id}/${sample.id}.libsize.txt"), emit: counts
    tuple val(sample), path("${stage}/${sample.id}/${sample.id}.is_paired.txt"), emit: is_paired
    
    script:
    """
    mkdir -p ${stage}/${sample.id}
    samtools flagstat $bam > "${stage}/${sample.id}/${sample.id}.flagstats.txt"
    head -n 1 "${stage}/${sample.id}/${sample.id}.flagstats.txt" | awk '{print \$1 + \$3}' > "${stage}/${sample.id}/${sample.id}.libsize.txt"
    grep -m 1 "paired in sequencing" "${stage}/${sample.id}/${sample.id}.flagstats.txt" | awk '{npaired = \$1 + \$3; if (npaired==0) {print "unpaired"} else {print "paired"};}' > "${stage}/${sample.id}/${sample.id}.is_paired.txt"
    """
}


process flagstats_libtype {
    container "quay.io/biocontainers/gawk:5.1.0--2"
    label "default"
    publishDir "${params.output_dir}", mode: "copy"

    input:
    path(files)

    output:
    path("stats/library_type.txt")

    script:
    """
    mkdir -p stats/
    find . -maxdepth 1 -mindepth 1 -name '*is_paired.txt' | xargs -I {} awk -v OFS='\t' '{ print gensub(/.+\\/(.+).is_paired.txt/, "\\\\1", "g", FILENAME), \$0;}' {} > stats/library_type.txt
    """
}


process collate_stats {
    // container "registry.git.embl.org/schudoma/portraits_metatraits:latest"
    container "quay.io/biocontainers/pandas:2.2.1"
    label "default"

    input:
    path(stats_files)

    output:
    path("reports/read_count_table.txt")

    script:
    """
    mkdir -p reports/
    collate_stats.py . > reports/read_count_table.txt
    """
}
