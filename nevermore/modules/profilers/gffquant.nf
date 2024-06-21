params.gq_aligner = "bwa"
params.gq_min_seqlen = params.min_alignment_length
params.gq_single_end_library = params.single_end_library
params.gq_min_identity = params.min_identity
params.gq_mode = "genes"
params.gq_ambig_mode = "1overN"


process stream_gffquant {
	publishDir params.output_dir, mode: "copy"
	label "gffquant"
	label "large"

	tag "gffquant.${sample}"

	input:
		tuple val(sample), path(fastqs)
		path(gq_db)
		// path(reference)
	output:
		tuple val(sample), path("profiles/${sample}/*.{txt.gz,pd.txt}"), emit: results //, optional: (!params.gq_panda) ? true : false
		tuple val(sample), path("profiles/${sample}/*.{txt.gz,pd.txt}"), emit: profiles //, optional: (params.gq_panda) ? true : false
		tuple val(sample), path("logs/${sample}.log")
		tuple val(sample), path("alignments/${sample}/${sample}*.sam"), emit: alignments, optional: true

	script:
			def gq_output = "-o profiles/${sample}/${sample}"

			def gq_params = "-m ${params.gq_mode} --ambig_mode ${params.gq_ambig_mode}"
			// gq_params += (params.gq_strand_specific) ? " --strand_specific" : ""
			gq_params += (params.gq_min_seqlen) ? (" --min_seqlen " + params.gq_min_seqlen) : ""
			gq_params += (params.gq_min_identity) ? (" --min_identity " + params.gq_min_identity) : ""
			// gq_params += (params.gq_restrict_metrics) ? " --restrict_metrics ${params.gq_restrict_metrics}" : ""
			// gq_params += (params.gq_keep_alignments) ? " --keep_alignment_file ${sample}.sam" : ""
			// gq_params += (params.gq_unmarked_orphans) ? " --unmarked_orphans" : ""

			gq_params += " -t ${task.cpus}"

			if (params.gq_mode == "domain") {
				gq_params += " --db_separator , --db_coordinates hmmer"
			}

			def input_files = ""
			// we cannot auto-detect SE vs. PE-orphan!
			if (params.gq_single_end_library) {
				//input_files += "--singles \$(find . -maxdepth 1 -type l -name '*_R1.fastq.gz')"	
				input_files += "--fastq-singles ${fastqs}"
			} else {
				r1_files = fastqs.findAll( { it.name.endsWith("_R1.fastq.gz") && !it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )
				r2_files = fastqs.findAll( { it.name.endsWith("_R2.fastq.gz") } )
				orphans = fastqs.findAll( { it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )

				if (r1_files.size() != 0) {
					input_files += "--fastq-r1 ${r1_files.join(' ')}"
				}
				if (r2_files.size() != 0) {
					input_files += " --fastq-r2 ${r2_files.join(' ')}"
				}
				if (orphans.size() != 0) {
					input_files += " --fastq-orphans ${orphans.join(' ')}"
				}

				// input_files += "--fastq-r1 \$(find . -maxdepth 1 -type l -name '*_R1.fastq.gz' | grep -v singles)"
				// input_files += " --fastq-r2 \$(find . -maxdepth 1 -type l -name '*_R2.fastq.gz')"
				// input_files += " --fastq-orphans \$(find . -maxdepth 1 -type l -name '*singles*.fastq.gz')"
			}
	
			def gq_cmd = "gffquant ${gq_output} ${gq_params} --db GQ_DATABASE --aligner ${params.gq_aligner} ${input_files}"
			def mkdir_alignments = (params.keep_alignment_file != null && params.keep_alignment_file != false) ? "mkdir -p alignments/${sample}/" : ""
			// --reference \$(readlink ${reference})
			// cp -v ${gq_db}/*sqlite3 GQ_DATABASE
			// ref=\$(ls ${gq_db}/*.bwt | sed "s/\.bwt//")
			"""
			set -e -o pipefail
			mkdir -p logs/ tmp/ profiles/
			${mkdir_alignments}
			echo 'Copying database...'
			cp -v \$(dirname \$(readlink ${gq_db}))/*sqlite3 GQ_DATABASE


			${gq_cmd} --reference \$(readlink ${gq_db}) &> logs/${sample}.log
			rm -rfv GQ_DATABASE* tmp/
			"""

}

process run_gffquant {
	publishDir params.output_dir, mode: "copy"
	label "gffquant"

	input:
	tuple val(sample), path(alignments) //, path(readcounts)
	path(gq_db)

	output:
	tuple val(sample), path("profiles/${sample}/*.txt.gz"), emit: results
	tuple val(sample), path("logs/${sample}.log")

	script:
	def gq_output = "-o profiles/${sample}/${sample}"

	def gq_params = "-m ${params.gq_mode} --ambig_mode ${params.gq_ambig_mode}"
	gq_params += (params.gq_strand_specific) ? " --strand_specific" : ""
	gq_params += (params.gq_unmarked_orphans) ? " --unmarked_orphans" : ""
	gq_params += (params.gq_min_seqlen) ? (" --min_seqlen " + params.gq_min_seqlen) : ""
	gq_params += (params.gq_min_identity) ? (" --min_identity " + params.gq_min_identity) : ""
	// gq_params += (params.bam_input_pattern) ? " --import_readcounts \$(grep -o '[0-9]\\+' ${readcounts})" : ""
	gq_params += (params.gq_restrict_metrics) ? " --restrict_metrics ${params.gq_restrict_metrics}" : ""
	// gq_params += (params.bam_input_pattern || !params.large_reference) ? (" --bam") : " --format sam"
	def formatted_input = (params.bam_input_pattern || !params.large_reference) ? "--bam ${alignments}" : "--sam ${alignments}"

	def gq_dbformat = (params.gq_mode == "domain") ? "--db_coordinates ${params.gq_db_coordinates} --db_separator ${params.gq_db_separator}" : ""
	def gq_cmd = "gffquant ${gq_output} ${gq_params} --db gq_db.sqlite3 ${gq_dbformat}"


	def mk_aln_sam = ""
	if (params.bam_input_pattern) {

		if (params.do_name_sort) {
			gq_cmd = "samtools collate -@ ${task.cpus} -O ${alignments} tmp/collated_bam | ${gq_cmd} --bam -"
		} else {
			gq_cmd = "${gq_cmd} ${formatted_input}"
		}

	} else if (params.large_reference) {

		mk_aln_sam += "echo 'Making alignment stream...'\n"
		if (alignments instanceof Collection && alignments.size() >= 2) {
			mk_aln_sam += "cat ${sample}.sam > tmp/alignments.sam \n"
			mk_aln_sam += "grep -v '^@' ${sample}.singles.sam >> tmp/alignments.sam"
		} else {
			mk_aln_sam += "ln -s ${alignments[0]} tmp/alignments.sam"
		}
		gq_cmd = "cat tmp/alignments.sam | ${gq_cmd} --sam -"

	} else {

		gq_cmd = "${gq_cmd} ${formatted_input}"

	}

	"""
	set -e -o pipefail
	mkdir -p logs/ tmp/ profiles/
	echo 'Copying database...'
	cp -v ${gq_db} gq_db.sqlite3
	${mk_aln_sam}
	${gq_cmd} &> logs/${sample}.log
	rm -rfv gq_db.sqlite3* tmp/
	"""
}

params.gq_collate_columns = "uniq_scaled,combined_scaled"

process collate_feature_counts {
	publishDir params.output_dir, mode: "copy"
	label "collate_profiles"
	label "gffquant"

	input:
	tuple val(sample), path(count_tables), val(column)
	val(suffix)

	output:
	path("collated/*.txt.gz"), emit: collated, optional: true

	script:
	def suffix_param = (suffix != "") ? "--suffix ${suffix}" : ""
	"""
	mkdir -p collated/

	collate_counts . -o collated/collated -c ${column} ${suffix_param}
	"""
}
