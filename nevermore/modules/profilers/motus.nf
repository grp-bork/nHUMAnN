params.motus_tax_level = "mOTU"
params.motus_min_length = 75
params.motus_n_marker_genes = 3


process motus {
    container "quay.io/biocontainers/motus:3.1.0--pyhdfd78af_0"

    input:
    tuple val(sample), path(fastqs)
	path(motus_db)

    output:
    tuple val(sample), path("${sample.id}/${sample.id}.motus.txt"), emit: motus_profile

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
	} else if (orphans.size() != 0) {
		input_files += "-s ${orphans.join(' ')}"
	}


    // def motus_input = (sample.is_paired) ? "-f ${sample.id}_R1.fastq.gz -r ${sample.id}_R2.fastq.gz" : "-s ${sample.id}_R1.fastq.gz";
    
    """
    mkdir -p ${sample.id}
    motus profile -n ${sample.id} -t $task.cpus -k ${params.motus_tax_level} -c -v 7 -q -l ${params.motus_min_length} -g ${params.motus_n_marker_genes} -db ${motus_db} ${input_files} > ${sample.id}/${sample.id}.motus.txt
    """
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