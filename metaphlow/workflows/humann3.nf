include { reduce_metaphlan_profiles; generate_humann_joint_index; run_humann3; humann_join_tables } from "../modules/profilers/humann3"


workflow humann3 {

	take:
		mp4_abundance_table
		mp4_tables
		fastqs

	main:
		reduce_metaphlan_profiles(
			mp4_abundance_table,
			"max"
		)

		generate_humann_joint_index(
			reduce_metaphlan_profiles.out.mp_reduced_profiles,
			params.humann_nuc_db
		)

		humann_input_ch = fastqs
			.join(mp4_tables, remainder: false)
			.map { sample, fastq, table -> return tuple(sample.id, table, fastq) }

		run_humann3(
			humann_input_ch,
			generate_humann_joint_index.out.joint_bt2_index,
			generate_humann_joint_index.out.chocophlan_db,
			params.humann_prot_db
		)

		collate_ch = run_humann3.out.hm_genefamilies
			.mix(run_humann3.out.hm_pathabundance)
			.mix(run_humann3.out.hm_pathcoverage)
			.mix(run_humann3.out.hm_genefamilies_relab)
			.mix(run_humann3.out.hm_pathabundance_relab)
			.mix(run_humann3.out.hm_table_stratified)
			.mix(run_humann3.out.hm_table_unstratified)
			.map { sample, file ->
				[file.name.replaceAll(/${sample}_/, "").replaceAll(/\.tsv$/, ""), file]
			}
			.groupTuple(by: 0, sort: true)

		humann_join_tables(collate_ch)

		// reformat_genefamily_table(run_humann3.out.hm_genefamilies) // now part of run_humann3
      
}


