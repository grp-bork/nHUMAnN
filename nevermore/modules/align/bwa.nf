process bwa_mem_align {
    container "registry.git.embl.de/schudoma/align-docker:latest"
    label 'align'

    input:
    tuple val(sample), path(reads)
    path(reference)
    val(do_name_sort)

    output:
    tuple val(sample), path("${sample.id}.bam"), emit: bam

    script:
    def maxmem = task.memory.toGiga()
    def align_cpus = 4 // figure out the groovy division garbage later (task.cpus >= 8) ?
    def sort_cpus = 4
    def blocksize = "-K 10000000"  // shamelessly taken from NGLess
    
    
    r1_files = reads.findAll( { it.name.endsWith("_R1.fastq.gz") && !it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )
    r2_files = reads.findAll( { it.name.endsWith("_R2.fastq.gz") } )
    orphan_files = reads.findAll( { it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )

    def r1_input = ""
    if (r1_files.size() != 0) {
        r1_input += "${r1_files.join(' ')}"
    } else if (orphan_files.size() != 0) {
        r1_input += "${orphan_files.join(' ')}"
    }
    def pre_sort_cmd_1 = "sortbyname.sh -Xmx${maxmem}g in=${r1_input} out=${sample.id}_R1.sorted.fastq.gz interleaved=f"
    def pre_sort_cmd_2 = ""
    def r2_input = ""
    def reads2 = ""
    if (r2_files.size() != 0) {
        r2_input += "${r2_files.join(' ')}"
        pre_sort_cmd_2 = "sortbyname.sh -Xmx${maxmem}g in=${r2_input} out=${sample.id}_R2.sorted.fastq.gz interleaved=f"
        reads2 = "${sample.id}_R2.sorted.fastq.gz"
    }

    def sort_cmd = (do_name_sort) ? "samtools collate -@ ${sort_cpus} -o ${sample.id}.bam - tmp/collated_bam" : "samtools sort -@ ${sort_cpus} -o ${sample.id}.bam -"

    def read_group_id = (sample.library == "paired") ? ((sample.is_paired) ? 2 : 2) : 1
    def read_group = "'@RG\\tID:${read_group_id}\\tSM:${sample.id}'"

    pre_sort_cmd_1 = ""
    pre_sort_cmd_2 = ""

    """
    set -e -o pipefail
    mkdir -p tmp/
    ${pre_sort_cmd_1}
    ${pre_sort_cmd_2}
    bwa mem -R ${read_group} -a -t ${align_cpus} ${blocksize} \$(readlink ${reference}) ${r1_input} ${r2_input} | samtools view -F 4 -buSh - | ${sort_cmd}
    rm -rvf tmp/ *.sorted.fastq.gz
    """
    // sortbyname.sh -Xmx${maxmem}g in=${sample.id}_R1.fastq.gz out=${sample.id}_R1.sorted.fastq.gz interleaved=f
}
