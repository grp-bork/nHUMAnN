// this nonsense parameter-rearrangement is necessary for nf-core schema/clowm compatibility

params.subsample_random_seed = 313
params.subsample = [:]
if (!params.subsample.random_seed) {
	params.subsample.random_seed = params.subsample_random_seed
}


process calculate_library_size_cutoff {
	label "tiny"


	input:
	path(readcounts)
	val(percentile)
	
	output:
	path("library_sizes.txt"), emit: library_sizes

	script:
	"""
	#!/usr/bin/env python3
	
	import glob
	import statistics

	"A_Niclo_T3h_r2.metaT.singles.raw.txt"
	d = dict(
		sorted(
			((f[:-8], int(open(f, "rt").read().strip().split("\t")[-1]))
			for f in glob.glob("*.raw.txt")),
			key=lambda x:x[1]
		)
	)
	percentile = ${percentile}
	try:
		percentiles = statistics.quantiles(d.values(), n=100)
	except statistics.StatisticsError:
		percentiles = None
	if percentiles is not None:
		mean_low_counts = statistics.mean(v for v in d.values() if v < percentiles[percentile - 1])
	else:
		mean_low_counts = list(d.values())[0]

	with open('library_sizes.txt', 'wt') as _out:
		print(*('sample', 'size', 'do_subsample', 'target_size'), sep='\\t', file=_out)
		for k, v in d.items():
			if percentiles is not None:
				do_subsample = v >= percentiles[percentile - 1]
			else:
				do_subsample = False
			print(k, v, int(do_subsample), int(mean_low_counts + 0.5), sep='\\t', file=_out)

	print(mean_low_counts)

	"""
	// nlibs=\$(cat ${readcounts} | wc -l)
	// cat ${readcounts} | sort -k1,1g | awk -v nlibs=\$nlibs 'BEGIN {q75=int(nlibs*0.75 + 0.5)} NR<q75 {print;}' | awk '{sum+=\$1} END {printf("%d\\n", sum/NR) }'
	// cat ${readcounts} | sort -k1,1g | awk -v nlibs=\$nlibs 'BEGIN {q75=int(nlibs*0.75 + 0.5)} NR<q75 {sum+=$1; n+=1;} END {printf("%d\n",sum/n) }'
}

process subsample_reads {
	container "quay.io/biocontainers/seqtk:1.4--he4a0461_2"
	label "medium"

	input:
	tuple val(sample), path(fastqs), val(target_size)

	output:
	tuple val(sample), path("subsampled/${sample.id}/*fastq.gz"), emit: subsampled_reads
	tuple val(sample), path("SUBSAMPLE.ok")

	script:

	def compression = (fastqs[0].name.endsWith(".gz")) ? "gz" : "bz2"
	def decomp = (compression == "gz") ? "gzip" : "bzip2"
	def seqtk_calls = ""
	def r1_files = fastqs.findAll( { it.name.endsWith("_R1.fastq.${compression}") } )
	def r2_files = fastqs.findAll( { it.name.endsWith("_R2.fastq.${compression}") } )
	    

	if (r1_files.size() != 0) {
        def r1_prefix = r1_files[0].name.replaceAll(/\.fastq.(gz|bz2)$/, "")
		seqtk_calls += "${decomp} -dc ${r1_files[0]} | seqtk sample -s ${params.subsample.random_seed} - ${target_size} | gzip -c - > subsampled/${sample.id}/${sample.id}_R1.fastq.gz\n"		
	}
	if (r2_files.size() != 0) {
        def r2_prefix = r2_files[0].name.replaceAll(/\.fastq.(gz|bz2)$/, "")
		seqtk_calls += "${decomp} -dc ${r2_files[0]} | seqtk sample -s ${params.subsample.random_seed} - ${target_size} | gzip -c - > subsampled/${sample.id}/${sample.id}_R2.fastq.gz\n"
	}

	"""
	set -e -o pipefail
	mkdir -p subsampled/${sample.id}/

	${seqtk_calls}

	touch SUBSAMPLE.ok
	"""


}

