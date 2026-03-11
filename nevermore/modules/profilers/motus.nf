params.motus_tax_level = "mOTU"
params.motus_min_length = 75
params.motus_n_marker_genes = 3
params.motus_readcount_type = "insert.scaled_counts"
params.motus_run_mapsnv = false
params.motus_full_rank_taxonomy = false
params.motus_print_counts = false

process run_motus {
    container "quay.io/biocontainers/motus:3.1.0--pyhdfd78af_0"

    input:
    tuple val(sample), path(fastqs)
	path(motus_db)

    output:
    tuple val(sample), path("${sample.id}/${sample.id}.motus.txt"), emit: motus_profile
    tuple val(sample), path("${sample.id}/${sample.id}.motus.bam"), emit: motus_bam

    script:

    def input_files = ""
	def r1_files = fastqs.findAll( { it.name.endsWith("_R1.fastq.gz") && !it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )
	def r2_files = fastqs.findAll( { it.name.endsWith("_R2.fastq.gz") } )
	def orphans = fastqs.findAll( { it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )

	if (r1_files.size() != 0 && r2_files.size() != 0) {
		input_files += "-f ${r1_files.join(' ')} -r ${r2_files.join(' ')}"
	} else if (r1_files.size() != 0) {
		input_files += "-s ${r1_files.join(' ')}"
	} else if (r2_files.size() != 0) {
		input_files += "-s ${r2_files.join(' ')}"
	}
    
    if (orphans.size() != 0) {
		input_files += " -s ${orphans.join(' ')}"
	}


    // def motus_input = (sample.is_paired) ? "-f ${sample.id}_R1.fastq.gz -r ${sample.id}_R2.fastq.gz" : "-s ${sample.id}_R1.fastq.gz";
    def mapsnv_cmd = ""
    if (params.motus_run_mapsnv) {
        mapsnv_cmd += "motus map_snv -t ${task.cpus} -db ${motus_db} ${input_files} > ${sample.id}/${sample.id}.motus.bam"
    } else {
        mapsnv_cmd += "touch ${sample.id}/${sample.id}.motus.bam"
    }

    def full_rank = ""
    if (params.motus_full_rank_taxonomy == true) {
        full_rank = "-q"
    }
    def print_counts = ""
    if (params.motus_print_counts) {
        print_counts = "-c"
    }
    
    """
    mkdir -p ${sample.id}
    motus profile -n ${sample.id} -t ${task.cpus} -k ${params.motus_tax_level} ${print_counts} -v 7 ${full_rank} -l ${params.motus_min_length} -g ${params.motus_n_marker_genes} -y ${params.motus_readcount_type} -db ${motus_db} ${input_files} > ${sample.id}/${sample.id}.motus.txt
    ${mapsnv_cmd}
    """

    // # generate profile
    // motus profile \
    // -f QC/${ID}.R1.fastq.gz -r QC/${ID}.R2.fastq.gz \
    // -g 1 \
    // -t 30 \
    // -y insert.raw_counts \
    // -o out_align/${ID}.profile.txt 
}



process motus_merge {
    container "quay.io/biocontainers/motus:3.1.0--pyhdfd78af_0"

    input:
    path(profiles)
    path(motus_db)

    output:
    path("motus_profiles/motus_merged.txt")

    script:
    """
    mkdir -p motus_profiles/ input/

    for f in ${profiles}; do ln -sf ../\$f input/; done

    motus merge -db ${motus_db} -d input/ -o motus_profiles/motus_merged.txt
    """


}