// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { SEQTK_SAMPLE              } from '../../../modules/nf-core/seqtk/sample/main'
include { FASTQ_ALIGN_DEDUP_BISMARK } from '../../../subworkflows/nf-core/fastq_align_dedup_bismark/main'
include { LAMBDACONVERSIONRATE      } from '../../../modules/local/lambda_conversion_rate/main'

workflow FASTQ_SUBSAMPLE_ALIGN_DEDUP_LAMBDA_CONVERSION_RATE_SEQTK_BISMARK {

    take:
    // TODO nf-core: edit input (take) channels
    reads         // [ meta, [ reads ] ]
    lambda_fasta  // path to lambda fasta file
    lambda_index  // path to lambda Bismark index
    sample_size   // number of reads to sample, default 1M

    main:
    // declare empty channels to be filled by module/workflow outputs
    ch_sample_reads           = Channel.empty()
    ch_bismark_summary        = Channel.empty()
    ch_lambda_conversion_rate = Channel.empty()
    ch_multiqc_files          = Channel.empty()
    ch_versions               = Channel.empty()


    // Subsample reads
    SEQTK_SAMPLE (
        reads,
        sample_size
    )
    ch_sample_reads    = SEQTK_SAMPLE.out.reads
    ch_versions        = ch_versions.mix(SEQTK_SAMPLE.out.versions)

    // Align subsampled reads to lambda genome and deduplicate
    FASTQ_ALIGN_DEDUP_BISMARK (
        ch_sample_reads,
        lambda_fasta,
        lambda_index,
        false,  // skip_deduplication = false
        false   // cytosine_report = false
    )
    ch_bismark_summary = FASTQ_ALIGN_DEDUP_BISMARK.out.bismark_summary
    ch_versions        = ch_versions.mix(FASTQ_ALIGN_DEDUP_BISMARK.out.versions)

    // Calculate lambda conversion rate
    LAMBDACONVERSIONRATE (
        ch_bismark_summary
    )
    ch_lambda_conversion_rate = LAMBDACONVERSIONRATE.out.conversion_rate
    ch_versions               = ch_versions.mix(LAMBDACONVERSIONRATE.out.versions)

    // Collect MultiQC inputs
    ch_multiqc_files = ch_lambda_conversion_rate.collect{ meta, rate -> rate }

    emit:
    // TODO nf-core: edit emitted channels
    lambda_conversion_rate = ch_lambda_conversion_rate // channel: [ val(meta), [ conversion_rate ] ]
    multiqc                = ch_multiqc_files          // path: *txt
    versions = ch_versions                             // path: version.txt
}

