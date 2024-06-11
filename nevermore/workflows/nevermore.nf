#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { nevermore_simple_preprocessing } from "./prep"
include { prepare_fastqs } from "../modules/converters/prepare_fastqs"
include { fastqc } from "../modules/qc/fastqc"
include { multiqc } from "../modules/qc/multiqc"
include { collate_stats } from "../modules/collate"
include { nevermore_align } from "./align"
include { nevermore_pack_reads } from "./pack"
include { nevermore_qa } from "./qa"
include { nevermore_decon } from "./decon"


def do_preprocessing = (!params.skip_preprocessing || params.run_preprocessing)
def do_alignment = params.run_gffquant || !params.skip_alignment
def do_stream = params.gq_stream


workflow nevermore_main {

	take:
		fastq_ch		

	main:
		if (do_preprocessing) {
	
			nevermore_simple_preprocessing(fastq_ch)
	
			preprocessed_ch = nevermore_simple_preprocessing.out.main_reads_out
			if (!params.drop_orphans) {
				preprocessed_ch = preprocessed_ch.mix(nevermore_simple_preprocessing.out.orphan_reads_out)
			}

			nevermore_decon(preprocessed_ch)
			preprocessed_ch = nevermore_decon.out.reads

		} else {
	
			preprocessed_ch = fastq_ch
	
		}
	
		nevermore_pack_reads(preprocessed_ch)

		collate_ch = Channel.empty()
		if (params.run_qa) {

			raw_counts_ch = (do_preprocessing) ? nevermore_simple_preprocessing.out.raw_counts : Channel.empty()

			nevermore_qa(
				nevermore_pack_reads.out.qa_fastqs,
				raw_counts_ch
			)

			collate_ch = nevermore_qa.out.readcounts_ch

		}

	emit:
		fastqs = nevermore_pack_reads.out.fastqs
		readcounts = collate_ch

}
