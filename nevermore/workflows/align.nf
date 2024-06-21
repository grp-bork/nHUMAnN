#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { bwa_mem_align } from "../modules/align/bwa"
include { minimap2_align } from "../modules/align/minimap2"
include { merge_and_sort; merge_sam } from "../modules/align/helpers"

def asset_dir = "${projectDir}/nevermore/assets"
def do_alignment = params.run_gffquant || !params.skip_alignment
def do_stream = params.gq_stream
def do_preprocessing = (!params.skip_preprocessing || params.run_preprocessing)

params.do_name_sort = true
params.align = [:]
params.align.run_minimap2 = false
params.align.run_bwa = false


workflow nevermore_align {

	take:
		fastq_ch

	main:

		alignment_ch = Channel.empty()
		aln_counts_ch = Channel.empty()

		/*	align merged single-read and paired-end sets against reference */

		if (params.align.run_minimap2) {
			minimap2_align(
				fastq_ch,
				params.reference,
				params.do_name_sort
			)

			minimap_aligned_ch = minimap2_align.out.sam
			.map { sample, sam ->
				sample_id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, "")
				return tuple(sample_id, sam)
			}
			.groupTuple(sort: true)

			/*	merge paired-end and single-read alignments into single per-sample bamfiles */
			merge_sam(minimap_aligned_ch
				.map { sample_id, samfiles ->
					def meta = [:]
					meta.id = sample_id
					return tuple(meta, samfiles)
			})

			alignment_ch = alignment_ch
				.mix(merge_sam.out.sam)			
			aln_counts_ch = aln_counts_ch
				.mix(merge_sam.out.flagstats)			
			
		}

		if (params.align.run_bwa) {
			bwa_mem_align(
				fastq_ch,
				params.reference,
				true
			)

			/*	merge paired-end and single-read alignments into single per-sample bamfiles */

			aligned_ch = bwa_mem_align.out.bam
				.map { sample, bam ->
					sample_id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, "")
					return tuple(sample_id, bam)
				}
				.groupTuple(sort: true)

			merge_and_sort(aligned_ch, true)

			alignment_ch = alignment_ch
				.mix(merge_and_sort.out.bam)
			aln_counts_ch = aln_counts_ch
				.mix(merge_and_sort.out.flagstats)
		}

	emit:
		alignments = alignment_ch 
		aln_counts = aln_counts_ch
}
