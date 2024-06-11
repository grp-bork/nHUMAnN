process sortmerna {
	container "quay.io/biocontainers/sortmerna:4.3.6--h9ee0642_0"

	input:
		tuple val(sample), path(fastqs)
		path(db)
	output:
		tuple val(sample), path("no_rrna/${sample.id}/*.fastq.gz"), emit: fastqs, optional: true
		tuple val(sample), path("rrna/${sample.id}/*.fastq.gz"), emit: rrna_fastqs, optional: true
	script:
		def mem_mb = task.memory.toMega()

		def reads = ""
		def mv_output = ""
		def pe_params = ""

		def r1_files = fastqs.findAll( { it.name.endsWith("_R1.fastq.gz") } )
		def r2_files = fastqs.findAll( { it.name.endsWith("_R2.fastq.gz") } )

		if (r1_files.size() != 0) {
			reads += "--reads ${r1_files[0]}"

			if (r2_files.size() != 0) {
				reads += " --reads ${r2_files[0]}"
				mv_output += "mv -v work/out/other_fwd.fq.gz no_rrna/${sample.id}/${sample.id}_R1.fastq.gz\n"
				mv_output += "mv -v work/out/other_rev.fq.gz no_rrna/${sample.id}/${sample.id}_R2.fastq.gz\n"
				pe_params += "--out2 --other --paired_in"

			} else {
				mv_output += "mv -v work/out/other.fq.gz no_rrna/${sample.id}/${sample.id}_R1.fastq.gz\n"	
			}
		}

	
		
		// if (sample.is_paired) {
		// 	reads += " --reads ${sample.id}_R2.fastq.gz"
		// 	// mv_output += "gzip -vc work/other_rev.fq > no_rrna/${sample.id}/${sample.id}_R2.fastq.gz\n"
		// 	mv_output += "mv -v work/out/other_fwd.fq.gz no_rrna/${sample.id}/${sample.id}_R1.fastq.gz\n"
		// 	mv_output += "mv -v work/out/other_rev.fq.gz no_rrna/${sample.id}/${sample.id}_R2.fastq.gz\n"
		// 	pe_params += "--out2 --other --paired_in"
		// } else {
		// 	mv_output = "mv -v work/out/other.fq.gz no_rrna/${sample.id}/${sample.id}_R1.fastq.gz"
		// }

		"""
		mkdir -p work/index/ no_rrna/${sample.id}/ rrna/${sample.id}/
		sortmerna --fastx --aligned --other --threads ${task.cpus} -m ${mem_mb} --workdir work/ --idx-dir \$(dirname \$(readlink ${db}))/index/ ${reads} ${pe_params} --ref ${db}

		${mv_output}
		"""

}



// def smk_get_sortmerna_params(wildcards, input):
//     """For snakemake: get SortMeRNA parameters"""
//     # reference files
//     # --fastx: output aligned reads
//     # -m: memory requirements for indexing
//     # --aligned/--other: output (non)-aligned reads
//     params = "{refs} --fastx -m {mem} --aligned --other".format(
//         refs=" ".join(["--ref %s" % ref for ref in input.refs]),
//         mem=int(MEMCORE[:-1]) * 1000
//     )
//     # additional parameters for PE reads
//     if hasattr(input, "r1"):
//         assert hasattr(input, "r2"), f"Found R1 but not R2 in {input}"
//         # --otu-map: OTU map (input to QIIME's make_otu_table.py)
//         # --paired_in: flags the PE reads as "aligned", when either of them is aligned
//         # --out2: output paired reads into separate files
//         params = f"{params} --otu_map --paired_in --out2"
//     return params

// rule filter_rrna_pe:
//     input:
//         r1="Preprocessing/mt.r1.trimmed.fq",
//         r2="Preprocessing/mt.r2.trimmed.fq",
//         refs=expand(
//             "{dbpath}/sortmerna/{fasta}.fasta",
//             dbpath=DBPATH, fasta=config["sortmerna"]["files"]
//         ),
//     output:
//         filt_r1="Preprocessing/mt.r1.trimmed.rna_filtered.fq",
//         filt_r2="Preprocessing/mt.r2.trimmed.rna_filtered.fq",
//         rrna_r1="Preprocessing/mt.r1.trimmed.rna.fq",
//         rrna_r2="Preprocessing/mt.r2.trimmed.rna.fq",
// 	    otu="Preprocessing/otu_map.txt",
//         # temp: working and index directories
//         wdir=temp(directory("Preprocessing/sortmerna_pe")),
//         idx=temp(directory("Preprocessing/sortmerna_idx")), # store index files separately so they can be re-used in filter_rrna_se
//     log:
//         "logs/preprocessing_filter_rrna_pe.log"
//     threads: getThreads(8)
//     resources:
//         runtime = "48:00:00",
//         mem = MEMCORE
//     params:
//         params=smk_get_sortmerna_params,
//     conda:
//         os.path.join(ENVDIR, "IMP_preprocessing.yaml")
//     message:
//         "filter_rna: filtering rRNA from mt PE reads."
//     shell:
//         """
//         if [ -d "{output.wdir}" ] && [ "$(ls -A {output.wdir})" ]; then
//             echo "ERROR: {output.wdir} exists and is not empty. Please remove it and re-run the pipeline."
//             exit 1
//         fi
        
//         """

// rule filter_rrna_se:
//     input:
//         se="Preprocessing/mt.se.trimmed.fq",
//         refs=expand(
//             "{dbpath}/sortmerna/{fasta}.fasta",
//             dbpath=DBPATH, fasta=config["sortmerna"]["files"]
//         ),
//         idx=rules.filter_rrna_pe.output.idx, # to avoid additional overhead
//     output:
//         filt_se="Preprocessing/mt.se.trimmed.rna_filtered.fq",
//         rrna_se="Preprocessing/mt.se.trimmed.rna.fq",
//         # temp: working directory
//         wdir=temp(directory("Preprocessing/sortmerna_se")),
    
//     params:
//         params=smk_get_sortmerna_params,
    
//     shell:
//         """
//         if [ -d "{output.wdir}" ] && [ "$(ls -A {output.wdir})" ]; then
//             echo "ERROR: {output.wdir} exists and is not empty. Please remove it and re-run the pipeline."
//             exit 1
//         fi
//         sortmerna {params.params} --threads {threads} --reads {input.se} --workdir {output.wdir}/ --idx-dir {input.idx}/ &>> {log}
//         cat {output.wdir}/out/aligned.log >> {log}
//         mv {output.wdir}/out/other.fq {output.filt_se}
//         mv {output.wdir}/out/aligned.fq {output.rrna_se}
//         """
