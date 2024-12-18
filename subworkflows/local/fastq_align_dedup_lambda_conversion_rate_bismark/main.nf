// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { SEQTK_SAMPLE } from '../../modules/nf-core/seqtk/sample/main'
include { SEQTK_SAMPLE } from '../../subworkflows/nf-core/fastq_align_dedup_bismark/main'
include { LAMBDACONVERSIONRATE } from '../../modules/local/lambda_conversion_rate/main'
include { SEQTK_SAMPLE } from '../../../modules/nf-core/seqtk/sample/main.nf'

workflow FASTQ_ALIGN_DEDUP_LAMBDA_CONVERSION_RATE_BISMARK {

    take:
    // TODO nf-core: edit input (take) channels
    reads         // [ meta, [ reads ] ]
    lambda_fasta  // path to lambda fasta file
    lambda_index  // path to lambda Bismark index
    subsample_size

    main:
    ch_subsample                  = Channel.empty()



    // Subsample reads
    SEQTK_SAMPLE (
        reads,
        subsample_size
    )
    ch_subsample = SEQTK_SAMPLE.out.fastq


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

