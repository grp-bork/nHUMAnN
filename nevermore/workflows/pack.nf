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
					sample.id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, ".singles")
					return tuple(sample, fastq)
			}
			.branch {
				single_end: it[0].library == "single"
				paired_end: it[0].library == "paired"
			}
		.set { single_reads_ch }

		def se_group_size = 2 - (params.drop_orphans ? 1 : 0)

		single_reads_ch.paired_end
			.groupTuple(sort: true, size: se_group_size, remainder: true)
			.branch {
				do_merge: it[1].size() > 1
				no_merge: true
			}
			.set { pe_singles_ch }

		merged_single_ch = pe_singles_ch.do_merge

		/*	then merge single-read file groups into single files */

		merge_single_fastqs(merged_single_ch)

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
			.mix(pe_singles_ch.no_merge)
			.mix(single_reads_ch.single_end)
			.mix(paired_ch)
			.mix(merge_single_fastqs.out.fastq)

		fastq_prep_ch = paired_ch
			.mix(single_reads_ch.single_end)
			.mix(pe_singles_ch.no_merge)
			.mix(merge_single_fastqs.out.fastq)

	emit:
		fastqs = fastq_prep_ch
		qa_fastqs = qa_ch
}
