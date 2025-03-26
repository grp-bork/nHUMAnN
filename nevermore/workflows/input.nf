nextflow.enable.dsl=2

include { classify_sample; classify_sample_with_library_info } from "../modules/functions"


params.bam_input_pattern = "**.bam"	

params.input_dir_structure = "tree"

def bam_suffix_pattern = params.bam_input_pattern.replaceAll(/\*/, "")


process transfer_fastqs {
	input:
		path(fastqs)
	output:
		path("fastq/*/**"), emit: fastqs
	script:
		"""
		transfer.py -i . -o fastq/
		"""
}


process transfer_bams {
	input:
		path(bamfiles)
	output:
		path("bam/*.bam"), emit: bamfiles
	script:
		"""
		mkdir -p bam/
		find . -maxdepth 1 -type l -name '*.bam' | xargs -I {} readlink {} | xargs -I {} rsync -vP {} bam/
		"""
}


process prepare_fastqs {
	label "default"

	input:
		tuple val(sample), path(files), val(remote_input), val(library_suffix)
	output:
		tuple val(sample), path("fastq/${sample}/${sample}*.fastq.{gz,bz2}"), emit: pairs, optional: true
		tuple val(sample), path("fastq/${sample}.singles/${sample}*.fastq.{gz,bz2}"), emit: singles, optional: true
		path("sample_library_info.txt"), emit: library_info

  script:
		def remote_option = (remote_input) ? "--remote-input" : ""
		def remove_suffix = (params.suffix_pattern) ? "--remove-suffix ${params.suffix_pattern}" : ""
		def input_dir_prefix = (params.input_dir) ? params.input_dir : params.remote_input_dir

		def custom_suffixes = (params.custom_fastq_file_suffixes) ? "--valid-fastq-suffixes ${params.custom_fastq_file_suffixes}" : ""

		def libsfx_param = (library_suffix != null) ? "--add_sample_suffix ${library_suffix}" : ""
		
		"""
		mkdir -p reads/${sample}
		ls ${files} | xargs -I{} sh -c 'ln -s ../../{} reads/${sample}/'
		prepare_fastqs.py -i reads -o fastq -p ${input_dir_prefix} ${custom_suffixes} ${remote_option} ${remove_suffix} ${libsfx_param}
		"""
}







workflow remote_fastq_input {
	take:
		fastq_ch

	main:

		transfer_fastqs(fastq_ch.collect())
		res_ch = transfer_fastqs.out.fastqs.flatten()

	emit:
		fastqs = res_ch
}


workflow remote_bam_input {
	take:
		bam_ch
	main:
		transfer_bams(bam_ch)
		res_ch = transfer_bams.out.bamfiles.flatten()
	emit:
		bamfiles = res_ch
}


workflow fastq_input {
	take:
		fastq_ch
		libsfx
	
	main:
		if (params.input_dir_structure == "flat") {
			fastq_ch = fastq_ch
				.map { file -> [ 
					file.getName()
						.replaceAll(/\.(fastq|fq)(\.(gz|bz2))?$/, "")
						.replaceAll(/[._]R?[12]$/, "")
						.replaceAll(/[._]singles$/, ""),
					file
				] }

		} else {
			fastq_ch = fastq_ch
				.map { file -> return tuple(file.getParent().getName(), file) }
		}

		fastq_ch = fastq_ch
			// .map { file -> return tuple(file.getParent().getName(), file) }
			.groupTuple(by: 0)
			.combine(libsfx)
			.map { sample_id, files, suffix -> return tuple(sample_id, files, (params.remote_input_dir != null || params.remote_input_dir), suffix) }

		if (params.ignore_samples) {
			ignore_samples = params.ignore_samples.split(",")
			print "Ignoring samples: ${ignore_samples}"
			fastq_ch = fastq_ch
				.filter { !ignore_samples.contains(it[0]) }
		}

		fastq_ch.dump(pretty: true, tag: "fastq_ch")
		prepare_fastqs(fastq_ch)
		prepare_fastqs.out.singles.mix(prepare_fastqs.out.pairs).dump(pretty: true, tag: "prepare_fastqs_out")

		library_info_ch = prepare_fastqs.out.library_info
			.splitCsv(header:false, sep:'\t', strip:true)
			.map { row ->
				return tuple(row[0], row[1])
			}

		prepped_fastq_ch = prepare_fastqs.out.singles
			.map { sample_id, files -> return tuple("${sample_id}.singles", files, false) }
			.mix(prepare_fastqs.out.pairs
				.map { sample_id, files -> return tuple(sample_id, files, true) }
			)
			.join(by: 0, library_info_ch)
			.map { sample_id, files, is_paired, library_is_paired ->
				def meta = [:]
				meta.id = sample_id
				meta.is_paired = is_paired
				meta.library = (library_is_paired == "1") ? "paired" : "single"
				return tuple(meta, [files].flatten())
			}
		prepped_fastq_ch.dump(pretty: true, tag: "prepped_fastq_ch")

	emit:
		fastqs = prepped_fastq_ch
}


workflow bam_input {
	take:
		bam_ch
	main:

		if (params.remote_input_dir) {
			bam_ch = remote_bam_input(bam_ch.collect())
		}

		bam_ch = bam_ch
			.map { file ->
				def sample = file.name.replaceAll(bam_suffix_pattern, "").replaceAll(/\.$/, "")
				return tuple(sample, file)
			}
			.groupTuple(sort: true)
			.map { classify_sample(it[0], it[1]) }

		fastq_ch = Channel.empty()
		if (params.do_bam2fq_conversion) {
			bam2fq(bam_ch)
			fastq_ch = bam2fq.out.reads
				.map { classify_sample(it[0].id, it[1]) }
		}
	emit:
		bamfiles = bam_ch
		fastqs = fastq_ch
}

