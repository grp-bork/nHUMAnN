process merge_and_sort {
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"
    label 'samtools'

    input:
    tuple val(sample), path(bamfiles)
    val(do_name_sort)

    output:
    tuple val(sample), path("bam/${sample.id}.bam"), emit: bam
    tuple val(sample), path("stats/bam/${sample.id}.flagstats.txt"), emit: flagstats

    script:
    def sort_order = (do_name_sort) ? "-n" : ""
    def merge_cmd = ""

    // need a better detection for this
    if (bamfiles instanceof Collection && bamfiles.size() >= 2) {
        merge_cmd = "samtools merge -@ $task.cpus ${sort_order} bam/${sample.id}.bam ${bamfiles}"
    } else {
        merge_cmd = "ln -s ../${bamfiles[0]} bam/${sample.id}.bam"
    }

    """
    mkdir -p bam/ stats/bam/
    ${merge_cmd}
    samtools flagstats bam/${sample.id}.bam > stats/bam/${sample.id}.flagstats.txt
    """
}


process merge_sam {
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"
    label 'samtools'

    input:
    tuple val(sample), path(samfiles)
    // val(do_name_sort)

    output:
    tuple val(sample), path("sam/${sample.id}.sam"), emit: sam
    tuple val(sample), path("stats/sam/${sample.id}.flagstats.txt"), emit: flagstats

    script:
    // def sort_order = (do_name_sort) ? "-n" : ""
    def merge_cmd = ""

    // need a better detection for this
    if (samfiles instanceof Collection && samfiles.size() >= 2) {
        // merge_cmd = "samtools merge -@ $task.cpus ${sort_order} bam/${sample.id}.bam ${bamfiles}"
        merge_cmd += "samtools view --no-PG -Sh ${samfiles[0]} > sam/${sample.id}.sam\n"
        merge_cmd += "samtools view -S ${samfiles[1]} >> sam/${sample.id}.sam"

    } else {
        merge_cmd = "ln -s ../${samfiles[0]} sam/${sample.id}.sam"
    }

    """
    mkdir -p sam/ stats/sam/
    ${merge_cmd}
    samtools flagstats sam/${sample.id}.sam > stats/sam/${sample.id}.flagstats.txt
    """
}


process db_filter {
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"
    label 'samtools'

    input:
    tuple val(sample), path(bam)
    path(db_bedfile)

    output:
    tuple val(sample), path("filtered_bam/${sample}.bam"), emit: bam
    tuple val(sample), path("stats/filtered_bam/${sample}.flagstats.txt"), emit: flagstats

    script:
    """
    mkdir -p filtered_bam/ stats/filtered_bam/
    bedtools intersect -u -ubam -a ${bam} -b ${db_bedfile} > filtered_bam/${sample}.bam
    samtools flagstats filtered_bam/${sample}.bam > stats/filtered_bam/${sample}.flagstats.txt
    """
}


process readcount {
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"
    label 'samtools'

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample}.readcounts.txt"), emit: readcounts

    script:
    """
    set -e -o pipefail
    mkdir -p tmp/
    samtools view ${sample}.bam | cut -f 1 | uniq | sort -u -T tmp/ | wc -l > ${sample}.readcounts.txt
    """
}


process db2bed3 {
    input:
    path(db)

    output:
    path("db.bed3"), emit: db

    script:
    """
    cp -v ${db} db.sqlite3
    sqlite3 db.sqlite3 'select seqid,start,end from annotatedsequence;' '.exit' | awk -F \\| -v OFS='\\t' '{print \$1,\$2-1,\$3}' | sort -k1,1 -k2,2g -k3,3g > db.bed3
    rm -v db.sqlite3*

    if [[ ! -s db.bed3 ]]; then
        echo 'Empty database!' >> /dev/stderr
        exit 1
    fi
    """
}
