#!/usr/bin/env nextflow

params.sample_id = null
params.fastq1 = null
params.fastq2 = null
params.benchmark_vcf = null
params.benchmark_bed = null
params.threads = 10
params.sample_type = "wgs" // wgs or wes
params.output_dir = null

// Fixed references
params.reference = "/mnt/disk1/namle/data/references/Homo_sapiens_assembly38.fasta"
params.known_sites = "/mnt/disk1/namle/data/references/Homo_sapiens_assembly38.known_indels.vcf.gz"

process germline {
    container "nvcr.io/nvidia/clara/clara-parabricks:4.5.0-1"
    containerOptions "--gpus '\"device=1,4\"'"
    publishDir "${params.output_dir}/${params.sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(fq1), path(fq2)
    path reference
    path known_sites

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam
    tuple val(sample_id), path("${sample_id}.germline.vcf"), emit: vcf
    path "${sample_id}.recal"

    script:
    """
    pbrun germline \\
        --ref ${reference} \\
        --x3 \\
        --low-memory \\
        --memory-limit 62 \\
        --htvc-low-memory \\
        --in-fq ${fq1} ${fq2} \\
        --knownSites ${known_sites} \\
        --out-bam ${sample_id}.bam \\
        --out-variants ${sample_id}.germline.vcf \\
        --out-recal-file ${sample_id}.recal
    """
}

process deepvariant {
    container "nvcr.io/nvidia/clara/clara-parabricks:4.5.0-1"
    containerOptions "--gpus '\"device=1,4\"'"
    publishDir "${params.output_dir}/${params.sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}.deepvariant.vcf")

    script:
    def wes_flag = params.sample_type == "wes" ? "--use-wes-model" : ""
    """
    pbrun deepvariant \\
        --ref ${reference} \\
        --in-bam ${bam} \\
        --num-streams-per-gpu 1 \\
        --out-variants ${sample_id}.deepvariant.vcf \\
        ${wes_flag}
    """
}

process benchmark {
    container "hap.py:latest"
    publishDir "${params.output_dir}/${params.sample_id}/benchmark", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf), val(prefix)
    path reference
    path reference_fai
    path benchmark_vcf
    path benchmark_vcf_tbi
    path benchmark_bed

    output:
    path "${sample_id}-${prefix}-*"

    script:
    """
    HGREF=${reference}
    /opt/hap.py/bin/hap.py \\
        ${benchmark_vcf} \\
        ${vcf} \\
        -f ${benchmark_bed} \\
        --threads ${params.threads} \\
        --engine vcfeval \\
        -r ${reference} \\
        -o ${sample_id}-${prefix}-vcfeval
    """
}

workflow {
    samples_ch = Channel.of([params.sample_id, file(params.fastq1), file(params.fastq2)])
    reference_ch = Channel.fromPath(params.reference).collect()
    reference_fai_ch = Channel.fromPath("${params.reference}.fai").collect()
    known_sites_ch = Channel.fromPath(params.known_sites).collect()
    benchmark_vcf_ch = Channel.fromPath(params.benchmark_vcf).collect()
    benchmark_vcf_tbi_ch = Channel.fromPath("${params.benchmark_vcf}.tbi").collect()
    benchmark_bed_ch = Channel.fromPath(params.benchmark_bed).collect()

    germline(samples_ch, reference_ch, known_sites_ch)
    deepvariant(germline.out.bam, reference_ch)

    // Combine VCF channels with prefixes for benchmarking
    vcf_ch = germline.out.vcf
        .map { sample_id, vcf -> [sample_id, vcf, "germline"] }
        .mix(deepvariant.out.map { sample_id, vcf -> [sample_id, vcf, "deepvariant"] })

    benchmark(vcf_ch, reference_ch, reference_fai_ch, benchmark_vcf_ch, benchmark_vcf_tbi_ch, benchmark_bed_ch)
}
