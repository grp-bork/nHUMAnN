process fastqc {
    container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"

    input:
    tuple val(sample), path(reads)
    val(stage)

    output:
    tuple val(sample), path("stats/${stage}/fastqc/*/*fastqc_data.txt"), emit: stats
    tuple val(sample), path("stats/${stage}/read_counts/*.${stage}.txt"), emit: counts

    script:
    def compression = (reads[0].name.endsWith(".gz")) ? "gz" : "bz2"
        
    def fastqc_cmd = "fastqc -t ${task.cpus} --extract --outdir=fastqc"
    def process_r2 = (sample.is_paired) ? "${fastqc_cmd} ${sample.id}_R2.fastq.${compression} && mv fastqc/${sample.id}_R2_fastqc/fastqc_data.txt fastqc/${sample.id}_R2_fastqc/${sample.id}_R2_fastqc_data.txt" : ""

    def fastqc_calls = ""
    def mv_calls = ""
	def r1_files = reads.findAll( { it.name.endsWith("_R1.fastq.${compression}") } ) // && !it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )
	def r2_files = reads.findAll( { it.name.endsWith("_R2.fastq.${compression}") } )
	// def orphans = fastqs.findAll( { it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )

    

	if (r1_files.size() != 0) {
        def r1_prefix = r1_files[0].name.replaceAll(/\.fastq.(gz|bz2)$/, "")
		fastqc_calls += "${fastqc_cmd} ${r1_files[0]}\n"
        mv_calls += "mv fastqc/${r1_prefix}_fastqc/fastqc_data.txt fastqc/${r1_prefix}_fastqc/${r1_prefix}_fastqc_data.txt\n"
	}
	if (r2_files.size() != 0) {
        def r2_prefix = r2_files[0].name.replaceAll(/\.fastq.(gz|bz2)$/, "")
		fastqc_calls += "${fastqc_cmd} ${r2_files[0]}\n"
        mv_calls += "mv fastqc/${r2_prefix}_fastqc/fastqc_data.txt fastqc/${r2_prefix}_fastqc/${r2_prefix}_fastqc_data.txt\n"
	}
	
    // ${fastqc_cmd} ${sample.id}_R1.fastq.${compression} && mv fastqc/${sample.id}_R1_fastqc/fastqc_data.txt fastqc/${sample.id}_R1_fastqc/${sample.id}_R1_fastqc_data.txt
    // ${process_r2}

    """
    set -e -o pipefail
    mkdir -p stats/${stage}/read_counts fastqc/

    ${fastqc_calls}
    ${mv_calls}

    grep "Total Sequences" fastqc/*/*data.txt > seqcount.txt
    echo \$(wc -l seqcount.txt)\$'\t'\$(head -n1 seqcount.txt | cut -f 2) > stats/${stage}/read_counts/${sample.id}.${stage}.txt
	mv fastqc stats/${stage}/
    """
}
