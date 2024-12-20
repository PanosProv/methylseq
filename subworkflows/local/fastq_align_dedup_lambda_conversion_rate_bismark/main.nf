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
    ch_sample_reads                  = Channel.empty()
    ch_versions                   = Channel.empty()



    // Subsample reads
    SEQTK_SAMPLE (
        reads,
        sample_size
    )
    ch_sample_reads    = SEQTK_SAMPLE.out.reads
    ch_versions = ch_versions.mix(SEQTK_SAMPLE.out.versions)

    // Align subsampled reads to lambda genome and deduplicate
    FASTQ_ALIGN_DEDUP_BISMARK (
        ch_sample_reads,
        lambda_fasta,
        lambda_index,
        false,  // skip_deduplication = false
        false   // cytosine_report = false
    )


    // TODO nf-core: substitute modules here for the modules of your subworkflow

    SAMTOOLS_SORT ( ch_bam )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    // TODO nf-core: edit emitted channels
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

