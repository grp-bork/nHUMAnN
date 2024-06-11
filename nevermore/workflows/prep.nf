#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { qc_bbduk } from "../modules/qc/bbduk"
include { qc_bbduk_stepwise_amplicon } from "../modules/qc/bbduk_amplicon"
include { qc_bbmerge } from "../modules/qc/bbmerge"
include { fastqc } from "../modules/qc/fastqc"
include { multiqc } from "../modules/qc/multiqc"
include { calculate_library_size_cutoff; subsample_reads } from "../modules/qc/subsample"


def merge_pairs = (params.merge_pairs || false)
def keep_orphans = (params.keep_orphans || false)

def asset_dir = "${projectDir}/nevermore/assets"

params.subsample = [:]

print asset_dir

process concat_singles {
    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}/${sample}.singles_R1.fastq.gz"), emit: reads

    script:
    """
    mkdir -p $sample
    cat ${reads} > ${sample}/${sample}.singles_R1.fastq.gz
    """
}


workflow nevermore_simple_preprocessing {

	take:

		fastq_ch

	main:
		rawcounts_ch = Channel.empty()
		if (params.run_qa || params.subsample.subset) {

			fastqc(fastq_ch, "raw")
			rawcounts_ch = fastqc.out.counts

			if (params.run_qa) {
				multiqc(
					fastqc.out.stats.map { sample, report -> report }.collect(),
					"${asset_dir}/multiqc.config",
					"raw"
				)
			}

			if (params.subsample.subset) {
				
				fastq_ch
					.branch {
						subsample: params.subsample.subset == "all" || it[0].library_source == params.subsample.subset
						no_subsample: true
					}
					.set { check_subsample_ch }
				// subsample_ch = fastq_ch
				// 	.filter { params.subsample.subset == "all" || it[0].library_source == params.subsample.subset }
				// subsample_ch.dump(pretty: true, tag: "subsample_ch")

				calculate_library_size_cutoff(
					fastqc.out.counts
						.filter { params.subsample.subset == "all" || it[0].library_source == params.subsample.subset }
						.map { sample, counts -> return counts }
						.collect(),
					params.subsample.percentile
				)
				calculate_library_size_cutoff.out.library_sizes.view()

				css_ch = check_subsample_ch.subsample
					.map { sample, fastqs -> return tuple(sample.id, sample, fastqs) }
					.join(
						
						calculate_library_size_cutoff.out.library_sizes
							.splitCsv(header: true, sep: '\t', strip: true)
							.map { row ->
								return tuple(row.sample, row.do_subsample == "1", row.target_size)
							},
							by: 0,
							remainder: true						
					)

				css_ch.dump(pretty: true, tag: "css_ch")


				// for some reason, .branch does not work here :S
				subsample_ch = css_ch
					.filter { it[3] }
					.map { sample_id, sample, fastqs, do_subsample, target_size ->
						return tuple(sample, fastqs, target_size)
					}
				subsample_ch.dump(pretty: true, tag: "subsample_ch")

				subsample_reads(subsample_ch)

				do_not_subsample_ch = css_ch
					.filter { !it[3] }
					.map { sample_id, sample, fastqs, do_subsample, target_size ->
						return tuple(sample, fastqs)
					}
					.mix(
						check_subsample_ch.no_subsample
					)
				do_not_subsample_ch.dump(pretty: true, tag: "do_not_subsample_ch")

				fastq_ch = do_not_subsample_ch
					.mix(subsample_reads.out.subsampled_reads)

				fastq_ch.dump(pretty: true, tag: "post_subsample_fastq_ch")
			}

		}

		processed_reads_ch = Channel.empty()
		orphans_ch = Channel.empty()

		if (params.amplicon_seq) {

			qc_bbduk_stepwise_amplicon(fastq_ch, "${asset_dir}/adapters.fa")
			processed_reads_ch = processed_reads_ch.mix(qc_bbduk_stepwise_amplicon.out.reads)
			orphans_ch = orphans_ch.mix(qc_bbduk_stepwise_amplicon.out.orphans)

		} else {

			qc_bbduk(fastq_ch, "${asset_dir}/adapters.fa")
			processed_reads_ch = processed_reads_ch.mix(qc_bbduk.out.reads)
			orphans_ch = qc_bbduk.out.orphans
				.map { sample, file -> 
					def meta = sample.clone()
					meta.id = sample.id + ".orphans"
					meta.is_paired = false
					return tuple(meta, file)
				}

		}

	emit:

		main_reads_out = processed_reads_ch
		orphan_reads_out = orphans_ch
		raw_counts = rawcounts_ch

}
