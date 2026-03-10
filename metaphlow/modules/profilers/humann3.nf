process reduce_metaphlan_profiles {
    container = "registry.git.embl.org/schudoma/humann3-docker:latest"
    label "humann3"
    label "medium"

	input:
		path(mp_collated_profiles)
		val(reduce_function)
	
	output:
		path("mp_${reduce_function}_reduced_profiles.txt"), emit: mp_reduced_profiles
	
	script:
		"""
		humann_reduce_table --input ${mp_collated_profiles} --output mp_${reduce_function}_reduced_profiles.txt.tmp --function ${reduce_function} --sort-by level
		sed -i "2 s/^/#/" mp_${reduce_function}_reduced_profiles.txt.tmp

		awk -v OFS='\\t' '{print \$1,\$2}' mp_${reduce_function}_reduced_profiles.txt.tmp > mp_${reduce_function}_reduced_profiles.txt
		"""
		// sed -i "2 s/^/xxx/" mp_${reduce_function}_reduced_profiles.txt.tmp
		// <(awk '{print \$1"\t"\$2"\tadditional_species"}' tmp) > joint_max_taxonomic_profile.tsv
		// mv mp_${reduce_function}_reduced_profiles.txt.tmp mp_${reduce_function}_reduced_profiles.txt		


}


process generate_humann_joint_index {
    container = "registry.git.embl.org/schudoma/humann3-docker:latest"
    label "humann3"
    label "process_high"
	
	input:
		path(mp_reduced_profiles)
		path(nuc_db)

	output:
		path("joint_bt2_index/**.{bt2,bt2l}"), emit: joint_bt2_index
		path("joint_bt2_index/**.ffn"), emit: chocophlan_db

	script:
		"""
		mkdir -p joint_bt2_index/
		echo -e "@dummy\nA\n+\n5" > joint.fastq

		humann --threads ${task.cpus} --input joint.fastq --input-format fastq --output joint_bt2_index/ --bypass-translated-search --taxonomic-profile ${mp_reduced_profiles} --nucleotide-database ${nuc_db}
		"""

}


process run_humann3 {
    publishDir params.output_dir, mode: "copy"
    container = "registry.git.embl.org/schudoma/humann3-docker:latest"
    label "humann3"
    label "process_high"

    input:
        tuple val(sample), path(mp_profile), path(fastq_files)
        path(joint_bt2_index)
		path(chocophlan_db)
		path(humann_prot_db)

    output:
        tuple val(sample), path("humann3/${sample}/${sample}_genefamilies.tsv"), emit: hm_genefamilies
        tuple val(sample), path("humann3/${sample}/${sample}_pathabundance.tsv"), emit: hm_pathabundance
        tuple val(sample), path("humann3/${sample}/${sample}_pathcoverage.tsv"), emit: hm_pathcoverage
        tuple val(sample), path("humann3/${sample}/${sample}_genefamilies.relab.tsv"), emit: hm_genefamilies_relab
        tuple val(sample), path("humann3/${sample}/${sample}_pathabundance.relab.tsv"), emit: hm_pathabundance_relab
        tuple val(sample), path("humann3/${sample}/${sample}_genefamilies.relab_stratified.tsv"), emit: hm_table_stratified
        tuple val(sample), path("humann3/${sample}/${sample}_genefamilies.relab_unstratified.tsv"), emit: hm_table_unstratified

    script:
    """
	mkdir -p humann3/${sample}/
    cat ${fastq_files} > merged.fq.gz

    humann \
    --taxonomic-profile ${mp_profile} \
    --nucleotide-database joint_bowtie2_index \
    --bypass-nucleotide-index \
    --protein-database ${humann_prot_db} \
    --input merged.fq.gz \
    --input-format fastq.gz \
    --output-basename ${sample} \
    --output humann3/${sample}/ \
    --threads ${task.cpus} \
    --remove-temp-output

    humann_renorm_table -i humann3/${sample}/${sample}_pathabundance.tsv -u relab -m community -s y -o humann3/${sample}/${sample}_pathabundance.relab.tsv
    humann_renorm_table -i humann3/${sample}/${sample}_genefamilies.tsv -u relab -m community -s y -o humann3/${sample}/${sample}_genefamilies.relab.tsv

    humann_split_stratified_table \
    --input humann3/${sample}/${sample}_genefamilies.relab.tsv \
    --output humann3/${sample}/

    rm merged.fq.gz
    """
}


process reformat_genefamily_table {
    publishDir params.output_dir, mode: "copy"
    label "humann3"
    label "process_single"
    container = "registry.git.embl.org/schudoma/humann3-docker:latest"

    input:
        tuple val(sample), path(hm_table)

    output:
        path "humann3/${sample}/${sample}_genefamilies.relab_stratified.tsv", emit: hm_table_stratified
        path "humann3/${sample}/${sample}_genefamilies.relab_unstratified.tsv", emit: hm_table_unstratified

    script:
    """
	mkdir -p humann3/${sample}/

    humann_renorm_table \
    --input ${hm_table} \
    --units relab \
    --output humann3/${sample}/${sample}_genefamilies.relab.tsv

    humann_split_stratified_table \
    --input ${sample}/${sample}_genefamilies.relab.tsv \
    --output ${sample}/
    """
}   

process humann_join_tables {
    tag "${tabletype}"
    publishDir params.output_dir, mode: "copy"
    container = "registry.git.embl.org/schudoma/humann3-docker:latest"
    label "humann3"
    label "process_single"
    
    input:
    tuple val(tabletype), path(tables)

    output:
    path("humann3/collated/humann3_${tabletype}.collated.tsv")

    script:
    """
    mkdir -p humann3/collated/

    humann_join_tables -i . -o humann3/collated/humann3_${tabletype}.collated.tsv --file_name .tsv
    """
    


}