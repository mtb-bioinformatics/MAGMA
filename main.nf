nextflow.enable.dsl = 2


//================================================================================
// Include sub-workflows/modules and (soft) override workflow-level parameters
//================================================================================

include { CALL_WF } from './workflows/call_wf.nf'
include { VALIDATE_FASTQS_WF } from './workflows/validate_fastqs_wf.nf'
include { MAP_WF } from './workflows/map_wf.nf'
include { MERGE_WF } from './workflows/merge_wf.nf'
include { MINOR_VARIANT_ANALYSIS_WF } from './workflows/minor_variant_analysis_wf.nf'
include { QUALITY_CHECK_WF } from './workflows/quality_check_wf.nf'
include { REPORTS_WF } from './workflows/reports_wf.nf'


//================================================================================
// Main workflow
//================================================================================

workflow {

    if (params.only_validate_fastqs) {

        VALIDATE_FASTQS_WF(params.input_samplesheet)

    } else {

        validated_reads_ch = VALIDATE_FASTQS_WF( params.input_samplesheet )

        QUALITY_CHECK_WF( validated_reads_ch )

        MAP_WF( validated_reads_ch )

        CALL_WF( MAP_WF.out.sorted_reads_ch )

        MINOR_VARIANT_ANALYSIS_WF(CALL_WF.out.reformatted_lofreq_vcfs_tuple_ch)

        MERGE_WF( CALL_WF.out.gvcf_ch,
                  CALL_WF.out.reformatted_lofreq_vcfs_tuple_ch, 
                  CALL_WF.out.cohort_stats_tsv,
                  MINOR_VARIANT_ANALYSIS_WF.out.approved_samples_ch,
                  MINOR_VARIANT_ANALYSIS_WF.out.rejected_samples_ch)

        REPORTS_WF(QUALITY_CHECK_WF.out.reports_fastqc_ch,
                   MINOR_VARIANT_ANALYSIS_WF.out.minor_variants_results_ch,
                   MERGE_WF.out.major_variants_results_ch)

    }

}

