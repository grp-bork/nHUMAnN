include { run_samestr_convert; run_samestr_merge; run_samestr_filter; run_samestr_stats; run_samestr_compare; run_samestr_summarize } from "../modules/profilers/samestr"
params.samestr_marker_db = "/scratch/schudoma/databases/samestr/mpa_vOct22_CHOCOPhlAnSGB_202212/marker_db/"

workflow samestr {

	take:
		mp4_sam
		mp4_tables

	main:
		run_samestr_convert(
			mp4_sam.join(mp4_tables),
			params.samestr_marker_db
		)

		grouped_npy_ch = run_samestr_convert.out.sstr_npy
			.join(run_samestr_convert.out.convert_sentinel, by: 0)
			.map { sample, data, sentinel -> return data }
			.flatten()
			.map { file ->
					def species = file.name.replaceAll(/[.].*/, "")
					return tuple(species, file)
			}
			.groupTuple(sort: true)
            
		run_samestr_merge(grouped_npy_ch, params.samestr_marker_db)
		run_samestr_filter(
			run_samestr_merge.out.sstr_npy,
			params.samestr_marker_db
		)
		run_samestr_stats(run_samestr_filter.out.sstr_npy, params.samestr_marker_db)
		run_samestr_compare(run_samestr_filter.out.sstr_npy, params.samestr_marker_db)

		run_samestr_summarize(
			run_samestr_compare.out.sstr_compare.collect(),
			mp4_tables.map { sample, table -> return table }.collect(),
			params.samestr_marker_db			
		)

}
