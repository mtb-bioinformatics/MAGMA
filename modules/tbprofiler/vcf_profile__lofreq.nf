process TBPROFILER_VCF_PROFILE__LOFREQ {
    tag "${params.vcf_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(name), path(mergedLofreqVcf)
        path(resistanceDb)

    output:
        tuple val(name), path("results/*")
        path("results/*"), emit: resistance_json


    script:
        def optionalDb  = resistanceDb ? "--db ${resistanceDb.name}" : ""

        """
        ${params.tbprofiler_path} vcf_profile \\
            ${optionalDb} \\
            ${mergedLofreqVcf}
        """

    stub:
        """
        mkdir results
        touch results/${sampleName}.results.json
        """

}
