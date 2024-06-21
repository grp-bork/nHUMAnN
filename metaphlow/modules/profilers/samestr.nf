process run_samestr_convert {
    container "registry.git.embl.de/schudoma/samestr-docker:latest"
    tag "${sample.id}"
    label "highmem_large"

    
    input:
		tuple val(sample), path(mp_sam), path(mp_profile)
		path(marker_db)

    output:
        tuple val(sample), path("sstr_convert/*/*.npz"), emit: sstr_npy, optional: true
        tuple val(sample), path("samestr_convert_DONE"), emit: convert_sentinel

    script:
    """
    set -e -o pipefail

    samestr --verbosity DEBUG \
    convert \
        --input-files ${mp_sam} \
        --min-vcov 1 \
        --min-aln-qual 0 \
        --marker-dir ${marker_db} \
        --output-dir sstr_convert/ \
        --nprocs ${task.cpus} \
        --tax-profiles-extension .txt

    touch samestr_convert_DONE
    """
}

process run_samestr_merge {
    publishDir params.output_dir, mode: "copy"
    container "registry.git.embl.de/schudoma/samestr-docker:latest"
    tag "${species}"
    label "highmem_large"
    
    input:
        tuple val(species), path(sstr_npy)
	path(marker_db)

    output:
        tuple \
            val(species), \
            path("sstr_merge/${species}.npz"), \
            path("sstr_merge/${species}.names.txt"), \
        emit: sstr_npy

    script:
    """
    samestr --verbosity DEBUG \
    merge \
        --input-files ${sstr_npy} \
        --output-dir sstr_merge/ \
	--marker-dir ${marker_db} \
        --clade ${species} \
        --nprocs ${task.cpus}
    """
}

process run_samestr_filter {
    container "registry.git.embl.de/schudoma/samestr-docker:latest"
    tag "${species}"
    label "highmem_large"
    
    input:
        tuple val(species), path(sstr_npy), path(sstr_names)
	path(marker_db)

    output:
        tuple \
            val(species), \
            path("sstr_filter/${species}.npz"), \
            path("sstr_filter/${species}.names.txt"), \
        emit: sstr_npy, optional: true

    script:
    // #    --global-pos-min-n-vcov 10 \
    // #    --sample-pos-min-n-vcov 2 \
    """

    samestr --verbosity DEBUG \
    filter \
        --input-files ${sstr_npy} \
        --input-names ${sstr_names} \
        --output-dir sstr_filter/ \
        --marker-dir ${marker_db} \
        --marker-trunc-len 50 \
        --global-pos-min-f-vcov 0.25 \
        --sample-pos-min-sd-vcov 3 \
        --samples-min-n-hcov 5000 \
        --sample-var-min-n-vcov 2 \
        --sample-var-min-f-vcov 0.025 \
        --clade-min-samples 1 \
        --nprocs ${task.cpus}
    """
}

process run_samestr_stats {
    publishDir params.output_dir, mode: "copy"
    container "registry.git.embl.de/schudoma/samestr-docker:latest"
    tag "${species}"
    label "large"
    
    input:
        tuple val(species), path(sstr_npy), path(sstr_names)
	path(marker_db)

    output:
        path "sstr_stats/${species}.aln_stats.txt", emit: sstr_stats

    script:
    """
    samestr --verbosity DEBUG \
    stats \
    --input-files ${sstr_npy} \
    --input-names ${sstr_names} \
    --marker-dir ${marker_db} \
    --nprocs ${task.cpus} \
    --output-dir sstr_stats/
    """
}

process run_samestr_compare {
    publishDir params.output_dir, mode: "copy"
    container "registry.git.embl.de/schudoma/samestr-docker:latest"
    tag "${species}"
    label "highmem_large"
    
    input:
        tuple val(species), path(sstr_npy), path(sstr_names)
	path(marker_db)

    output:
        tuple \
            path("sstr_compare/${species}.closest.txt"), \
            path("sstr_compare/${species}.fraction.txt"), \
            path("sstr_compare/${species}.overlap.txt"), \
        emit: sstr_compare

    script:
    """
    samestr --verbosity DEBUG \
    compare \
        --input-files ${sstr_npy} \
        --input-names ${sstr_names} \
        --marker-dir ${marker_db} \
        --output-dir sstr_compare/ \
        --nprocs ${task.cpus}
    """
}

process run_samestr_summarize {
    publishDir params.output_dir, mode: "copy"
    container "registry.git.embl.de/schudoma/samestr-docker:latest"
    label "large"
    
    input:
        path(sstr_data)
        path(mp_profiles)
	path(marker_db)

    output:
        tuple \
            path("sstr_summarize/taxon_counts.tsv"), \
            path("sstr_summarize/sstr_cooccurrences.tsv"), \
            path("sstr_summarize/sstr_strain_events.tsv"), \
        emit: sstr_summarize

    script:
    """
    mkdir profiles/
    # TODO: this will only work with metaphlan, not motus
    find . -maxdepth 1 -name '*.mp4.txt' -exec mv -v {} profiles/ \\;

    samestr --verbosity DEBUG \
    summarize \
        --input-dir ./ \
        --marker-dir ${marker_db} \
        --tax-profiles-dir ./profiles/ \
        --tax-profiles-extension .txt \
        --output-dir sstr_summarize/
    """
}
