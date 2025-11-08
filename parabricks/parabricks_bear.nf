#!/usr/bin/env nextflow

params.sample_id = null
params.fastq1 = null
params.fastq2 = null
params.benchmark_vcf = null
params.benchmark_bed = null
params.threads = 20
params.sample_type = "wgs" // wgs or wes
params.output_dir = null

// Fixed references
params.reference = "/mnt/nasdev2/namle/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
params.known_sites = "/mnt/nasdev2/namle/references/Homo_sapiens_assembly38.known_indels.vcf.gz"

process germline {
    container "nvcr.io/nvidia/clara/clara-parabricks:4.5.0-1"
    containerOptions "--gpus '\"device=3,4\"'"
    publishDir "${params.output_dir}/${params.sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(fq1), path(fq2)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam
    tuple val(sample_id), path("${sample_id}.germline.vcf"), emit: vcf
    path "${sample_id}.recal"

    script:
    """
    pbrun germline \\
        --ref ${params.reference} \\
        --in-fq ${fq1} ${fq2} \\
        --knownSites ${params.known_sites} \\
        --out-bam ${sample_id}.bam \\
        --out-variants ${sample_id}.germline.vcf \\
        --out-recal-file ${sample_id}.recal
    """
}

process deepvariant {
    container "nvcr.io/nvidia/clara/clara-parabricks:4.5.0-1"
    containerOptions "--gpus '\"device=3,4\"'"
    publishDir "${params.output_dir}/${params.sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.deepvariant.vcf")

    script:
    def wes_flag = params.sample_type == "wes" ? "--use-wes-model" : ""
    """
    pbrun deepvariant \\
        --ref ${params.reference} \\
        --in-bam ${bam} \\
        --out-variants ${sample_id}.deepvariant.vcf \\
        ${wes_flag}
    """
}

process benchmark {
    container "hap.py:latest"
    publishDir "${params.output_dir}/${params.sample_id}/benchmark", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf), val(prefix)

    output:
    path "${sample_id}-${prefix}-*"

    script:
    """
    /opt/hap.py/bin/hap.py \\
        ${params.benchmark_vcf} \\
        ${vcf} \\
        -f ${params.benchmark_bed} \\
        --threads ${params.threads} \\
        --engine vcfeval \\
        -r ${params.reference} \\
        -o ${sample_id}-${prefix}-vcfeval
    """
}

workflow {
    samples_ch = Channel.of([params.sample_id, file(params.fastq1), file(params.fastq2)])

    germline(samples_ch)
    deepvariant(germline.out.bam)

    // Combine VCF channels with prefixes for benchmarking
    vcf_ch = germline.out.vcf
        .map { sample_id, vcf -> [sample_id, vcf, "germline"] }
        .mix(deepvariant.out.map { sample_id, vcf -> [sample_id, vcf, "deepvariant"] })

    benchmark(vcf_ch)
}
