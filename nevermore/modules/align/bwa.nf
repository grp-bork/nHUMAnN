process bwa_mem_align {
    cpus 4
    memory { 20.GB * task.attempt }
    time { 3.d * task.attempt }
    container "registry.git.embl.org/schudoma/align-docker:latest"
    label 'align'

    input:
    tuple val(sample), path(reads)
    path(reference)
    val(do_name_sort)

    output:
    tuple val(sample), path("${sample.id}.bam"), emit: bam

    script:
    def maxmem = task.memory.toGiga()
    
    def align_cpus = 1
    def sort_cpus = 1
    if (task.cpus > 1) {
        def half_cpus = task.cpus.intdiv(2)
        sort_cpus = half.cpus
        align_cpus = task.cpus - half.cpus
    }
    def blocksize = "-K 10000000"  // shamelessly taken from NGLess
    
    // def sort_cmd = "samtools collate -@ ${sort_cpus} -o ${sample.id}.bam - tmp/collated_bam"
    
    r1_files = reads.findAll( { it.name.endsWith("_R1.fastq.gz") && !it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )
    r2_files = reads.findAll( { it.name.endsWith("_R2.fastq.gz") } )
    orphan_files = reads.findAll( { it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )

    def r1_input = ""
    if (r1_files.size() != 0) {
        r1_input += "${r1_files.join(' ')}"
    } else if (orphan_files.size() != 0) {
        r1_input += "${orphan_files.join(' ')}"
    }
    def r2_input = ""
    if (r2_files.size() != 0) {
        r2_input += "${r2_files.join(' ')}"
    }

    def sort_cmd = (do_name_sort) ? "samtools collate -@ ${sort_cpus} -o ${sample.id}.bam - tmp/collated_bam" : "samtools sort -@ ${sort_cpus} -o ${sample.id}.bam -"

    def read_group_id = (sample.library == "paired") ? ((sample.is_paired) ? 2 : 2) : 1
    def read_group = "'@RG\\tID:${read_group_id}\\tSM:${sample.id}'"

    """
    set -e -o pipefail
    mkdir -p tmp/
    bwa mem -R ${read_group} -a -t ${align_cpus} ${blocksize} \$(readlink ${reference}) ${r1_input} ${r2_input} | samtools view -buSh - | ${sort_cmd}
    rm -rvf tmp/ *.sorted.fastq.gz
    """
}
