include { reduce_metaphlan_profiles; generate_humann_joint_index; run_humann3 } from "../modules/profilers/humann3"


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

		// reformat_genefamily_table(run_humann3.out.hm_genefamilies) // now part of run_humann3
      
}


