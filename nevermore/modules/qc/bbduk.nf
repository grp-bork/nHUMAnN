params.qc_params_shotgun = "qtrim=rl trimq=3 maq=25 ktrim=r k=23 mink=11 hdist=1 ftm=5 entropy=0.5 entropywindow=50 entropyk=5 tpe tbo"
params.qc_minlen = 45

process qc_bbduk {
    container "quay.io/biocontainers/bbmap:39.06--h92535d8_0"
	label 'bbduk'
    tag "${sample.id}"


    input:
    tuple val(sample), path(reads)
	path(adapters)

    output:
    tuple val(sample), path("qc_reads/${sample.id}/${sample.id}_R*.fastq.gz"), emit: reads
    tuple val(sample), path("qc_reads/${sample.id}/${sample.id}.orphans_R1.fastq.gz"), emit: orphans, optional: true
    path("stats/qc/bbduk/${sample.id}.bbduk_stats.txt")
    tuple val(sample), path("qc_reads/${sample.id}/BBDUK_FINISHED"), emit: sentinel

    script:
    // def maxmem = task.memory.toGiga().intdiv(2) 
    def maxmem = task.memory.toGiga() 
    def compression = (reads[0].name.endsWith("gz")) ? "gz" : "bz2"

    def read2 = ""
    def orphan_check = ""

    def bb_params = params.qc_params_shotgun //.replaceAll(/maq=([0-9]+)/, "")
    
    def trim_params = "${bb_params} ref=${adapters} minlen=${params.qc_minlen}"

    def orphan_filter = ""

    def r1_files = reads.findAll( { it.name.endsWith("_R1.fastq.${compression}") } )
	def r2_files = reads.findAll( { it.name.endsWith("_R2.fastq.${compression}") } )

    def qenc_str = (params.phred64 != null && params.phred64 != false) ? "qin=64" : ""

    def read1 = ""
    def orphans = ""
    if (r1_files.size() != 0) {
        read1 += "in1=${r1_files[0]} out1=qc_reads/${sample.id}/${sample.id}_R1.fastq.gz"
        // read1 = "in1=${sample.id}_R1.fastq.${compression} out1=qc_reads/${sample.id}/${sample.id}_R1.fastq.gz"
        if (r2_files.size() != 0) {
            read2 += "in2=${r2_files[0]} out2=qc_reads/${sample.id}/${sample.id}_R2.fastq.gz outs=tmp_orphans.fq"
            orphans += "qc_reads/${sample.id}/${sample.id}.orphans_R1.fastq.gz"
            orphan_filter += "bbduk.sh -Xmx${maxmem}g t=${task.cpus} ${trim_params} in=tmp_orphans.fq out=${orphans}"
            orphan_check = """
            if [[ -z "\$(gzip -dc ${orphans} | head -n 1)" ]]; then
                rm ${orphans}
            fi
            """
        }
    }
	

    
    // if (sample.is_paired) {
    //     def orphans = "qc_reads/${sample.id}/${sample.id}.orphans_R1.fastq.gz"
    //     // read2 = "in2=${sample.id}_R2.fastq.${compression} out2=qc_reads/${sample.id}/${sample.id}_R2.fastq.gz outs=tmp_orphans.fq"
    //     // orphan_filter = "bbduk.sh -Xmx${maxmem}g t=${task.cpus} ${trim_params} in=tmp_orphans.fq out=${orphans}"

    //     orphan_check = """
    //     if [[ -z "\$(gzip -dc ${orphans} | head -n 1)" ]]; then
	// 		rm ${orphans}
	// 	fi
    //     """
    // }

    // def read1 = "in1=${sample.id}_R1.fastq.${compression} out1=qc_reads/${sample.id}/${sample.id}_R1.fastq.gz"
    
    def stats_out = "stats=stats/qc/bbduk/${sample.id}.bbduk_stats.txt"

    """
    set -e -o pipefail

    mkdir -p qc_reads/${sample.id}/ stats/qc/bbduk/
    bbduk.sh -Xmx${maxmem}g t=${task.cpus} ${trim_params} ${qenc_str} ${stats_out} ${read1} ${read2}
    ${orphan_filter}
    ${orphan_check}

    touch qc_reads/${sample.id}/BBDUK_FINISHED
    rm -vf *.fq
    """
}
