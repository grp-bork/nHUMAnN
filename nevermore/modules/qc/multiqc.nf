process multiqc {
    container "quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0"

    input:
    path(reports)
	path(multiqc_config)
	val(stage)

    output:
    path("reports/${stage}.multiqc_report.html")

    script:
    """
	mkdir -p reports/
    multiqc -o reports/ -n ${stage}.multiqc_report.html -c ${multiqc_config} .
    """
}
