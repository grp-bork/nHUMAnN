include { merge_single_fastqs } from "../modules/converters/merge_fastqs"


workflow nevermore_pack_reads {

	take:
		fastq_ch
	
	main:

		/* re-add pair information, which might have been lost upstream */

		fastq_ch = fastq_ch
			.map { sample, fastqs ->
			def meta = sample.clone()
			meta.is_paired = [fastqs].flatten().size() == 2
			return tuple(meta, fastqs)
		}

		fastq_ch.dump(pretty: true, tag: "pack_fastq_ch")

		/*	route all single-read files into a common channel */

		single_ch = fastq_ch
			.filter { it[0].is_paired == false }
			.map { sample, fastq ->
				def meta = sample.clone()
				meta.id = fastq.name.replaceAll(/_R1.fastq.gz$/, "")
				meta.merged = false
				return tuple(meta, fastq)
			}

		/*	route all paired-end read files into a common channel */

		paired_ch = fastq_ch
			.filter { it[0].is_paired == true }
			.map { sample, fastq ->
				def meta = sample.clone()
				meta.merged = true
				return tuple(meta, fastq)
			}

		/*	group all single-read files by sample and route into merge-channel */

		single_ch
			.map {
				sample, fastq ->
					def meta = sample.clone()
					meta.id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, ".singles")
					return tuple(meta, fastq)
			}
			.branch {
				single_end: it[0].library == "single"
				paired_end: it[0].library == "paired"
			}
		.set { single_reads_ch }

		def orphan_merge = !params.single_end_libraries && !params.drop_orphans && params.run_preprocessing //&& params.remove_host;
		def se_group_size = 2 - ((orphan_merge) ? 0 : 1);

		single_reads_ch.paired_end
			.branch {
				do_merge: it[0].multilib
				no_merge: true
			}
			.set { pe_singles_ch }

		merged_single_ch = pe_singles_ch.do_merge
			.map { meta, fastq -> [ meta.id, fastq ] }
			.groupTuple(by: 0, sort: true, size: se_group_size, remainder: true)
			.map { sample_id, fastqs ->
				def meta = [:]
				meta.id = sample_id
				meta.is_paired = false
				meta.library = "paired"
				meta.merged = true
				meta.multilib = true

				return [meta, [fastqs].flatten()]
			}

		merged_single_ch.dump(pretty: true, tag: "merged_single_ch")

		/*	then merge single-read file groups into single files */

		merged_ch = Channel.empty()
		if (!params.single_end_libraries) {

			merge_single_fastqs(merged_single_ch)
			merged_ch = merge_single_fastqs.out.fastq

		} 

		/* 	take all single-read files except for the qc-survivors,
			concat with merged single-read files (takes care of single-end qc-survivors),
			concat with paired-end files,
			and route them into a channel for post-qc fastqc analysis
			(THIS IS JUST FOR QA/STATS PURPOSES, NO WORRIES, THE SE READS ARE PROCESSED PROPERLY!)
		*/

		qa_ch = single_ch
			.filter { ! it[0].id.endsWith(".singles") }
			.map { sample, fastq ->
				def meta = sample.clone()
				meta.id = fastq.name.replaceAll(/_R1.fastq.gz$/, "")
				meta.merged = false
				return tuple(meta, fastq)
			}
			.mix(pe_singles_ch.no_merge)         // raw PE library orphans
			.mix(single_reads_ch.single_end)     // SE library reads
			.mix(paired_ch)                      // PE library pairs
			.mix(merged_ch)                      // merged preprocessed PE library orphans

		fastq_prep_ch = paired_ch
			.mix(single_reads_ch.single_end)
			.mix(pe_singles_ch.no_merge)
			.mix(merged_ch)

	emit:
		fastqs = fastq_prep_ch
		qa_fastqs = qa_ch
}
