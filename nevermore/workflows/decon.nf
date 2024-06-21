include { remove_host_kraken2_individual; remove_host_kraken2 } from "../modules/decon/kraken2"
include { hostile } from "../modules/decon/hostile"
include { sortmerna } from "../modules/decon/sortmerna"


workflow nevermore_decon {

		take:
			fastq_ch
		
		main:
			preprocessed_ch = Channel.empty()

			if (params.run_sortmerna) {

				fastq_ch
					.branch {
						metaT: it[0].containsKey("library_source") && it[0].library_source == "metaT"
						metaG: true
					}
					.set { for_sortmerna_ch }

				sortmerna(for_sortmerna_ch.metaT, params.sortmerna_db)
				preprocessed_ch = for_sortmerna_ch.metaG
					.mix(sortmerna.out.fastqs)

			} else {

				preprocessed_ch = fastq_ch

			}
	
			if (params.remove_host == "hostile") {

				hostile(preprocessed_ch, params.hostile_db)
				preprocessed_ch = hostile.out.reads

			} else if ((params.remove_host != false && params.remove_host != null ) || params.remove_host == "kraken") {

				remove_host_kraken2_individual(preprocessed_ch, params.remove_host_kraken2_db)	
				preprocessed_ch = remove_host_kraken2_individual.out.reads

			}

		emit:
			reads = preprocessed_ch

}
