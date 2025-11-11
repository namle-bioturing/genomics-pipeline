#!/usr/bin/env nextflow

params.sample_id = null
params.input = null  // Comma-separated: for fastq: R1_1,R2_1,R1_2,R2_2,... or for bam: file.bam,file.bam.bai
params.input_type = "fastq"  // fastq or bam
params.target_bed = null
params.benchmark_vcf = null
params.benchmark_bed = null
params.threads = 20
params.sample_type = "wgs" // wgs or wes
params.pb_deepvariant_model = null
params.gg_deepvariant_model = null
params.output_dir = null
params.interval = "-L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY -L chrM"

// Fixed references
params.reference = "/mnt/disk1/namle/data/references/Homo_sapiens_assembly38.fasta"
params.known_sites = "/mnt/disk1/namle/data/references/Homo_sapiens_assembly38.known_indels.vcf.gz"

process bam2fq {
    container "nvcr.io/nvidia/clara/clara-parabricks:4.5.0-1"
    containerOptions "--gpus '\"device=1,4\"'"
    publishDir "${params.output_dir}/${params.sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}_1.fastq.gz"), path("${sample_id}_2.fastq.gz")

    script:
    """
    pbrun bam2fq \\
        --in-bam ${bam} \\
        --num-threads ${params.threads} \\
        --out-prefix ${sample_id}
    """
}

process bamsort {
    container "nvcr.io/nvidia/clara/clara-parabricks:4.5.0-1"
    containerOptions "--gpus '\"device=1,4\"'"
    publishDir "${params.output_dir}/${params.sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam")

    script:
    """
    pbrun bamsort \\
        --x3 \\
        --gpusort \\
        --gpuwrite \\
        --ref ${reference} \\
        --num-zip-threads 16 \\
        --num-sort-threads 10 \\
        --mem-limit 62 \\
        --in-bam ${bam} \\
        --out-bam ${sample_id}.sorted.bam \\
        --sort-order queryname
    """
}

process markdup {
    container "nvcr.io/nvidia/clara/clara-parabricks:4.5.0-1"
    containerOptions "--gpus '\"device=1,4\"'"
    publishDir "${params.output_dir}/${params.sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}.markdup.bam"), path("${sample_id}.markdup.bam.bai")

    script:
    """
    pbrun markdup \\
        --x3 \\
        --ref ${reference} \\
        --mem-limit 62 \\
        --in-bam ${bam} \\
        --out-bam ${sample_id}.markdup.bam
    """
}

process bqsr {
    container "nvcr.io/nvidia/clara/clara-parabricks:4.5.0-1"
    containerOptions "--gpus '\"device=1,4\"'"
    publishDir "${params.output_dir}/${params.sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    path known_sites

    output:
    tuple val(sample_id), path(bam), path("${sample_id}.recal")

    script:
    """
    pbrun bqsr \\
        --ref ${reference} \\
        --in-bam ${bam} \\
        --knownSites ${known_sites} \\
        --out-recal-file ${sample_id}.recal
    """
}

process haplotypecaller {
    container "nvcr.io/nvidia/clara/clara-parabricks:4.5.0-1"
    containerOptions "--gpus '\"device=1,4\"'"
    publishDir "${params.output_dir}/${params.sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    path reference
    path recal_file

    output:
    tuple val(sample_id), path("${sample_id}.haplotypecaller.vcf")

    script:
    """
    pbrun haplotypecaller \\
        --x3 \\
        --htvc-low-memory \\
        --ref ${reference} \\
        --in-bam ${bam} \\
        --in-recal-file ${recal_file} \\
        ${params.interval} \\
        --out-variants ${sample_id}.haplotypecaller.vcf
    """
}

process germline {
    container "nvcr.io/nvidia/clara/clara-parabricks:4.5.0-1"
    containerOptions "--gpus '\"device=1,4\"'"
    publishDir "${params.output_dir}/${params.sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(fq_files)
    path reference
    path known_sites

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: bam
    tuple val(sample_id), path("${sample_id}.germline.vcf"), emit: vcf
    path "${sample_id}.recal"

    script:
    def fq_pairs = []
    if (fq_files instanceof List) {
        // Group files into pairs (R1, R2) and create separate --in-fq for each pair
        for (int i = 0; i < fq_files.size(); i += 2) {
            fq_pairs.add("--in-fq ${fq_files[i]} ${fq_files[i+1]}")
        }
    } else {
        fq_pairs.add("--in-fq ${fq_files}")
    }
    def in_fq_args = fq_pairs.join(' \\\n')
    """
    pbrun germline \
        --ref ${reference} \
        --x3 \
        --low-memory \
        --memory-limit 62 \
        --htvc-low-memory \
        ${in_fq_args} \
        --knownSites ${known_sites} \
        --out-bam ${sample_id}.bam \
        --out-variants ${sample_id}.germline.vcf \
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
    path model_file

    output:
    tuple val(sample_id), path("${sample_id}.deepvariant.vcf")

    script:
    def wes_flag = params.sample_type == "wes" ? "--use-wes-model" : ""
    def model_flag = model_file.name != 'NO_MODEL_FILE' ? "--bp-model-file ${model_file}" : ""
    """
    pbrun deepvariant \\
        --ref ${reference} \\
        --in-bam ${bam} \\
        --num-streams-per-gpu 1 \\
        --out-variants ${sample_id}.deepvariant.vcf \\
        ${params.interval} \\
        ${wes_flag} \\
        ${model_flag}
    """
}

process google_deepvariant {
    container "google/deepvariant:1.9.0"
    publishDir "${params.output_dir}/${params.sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    path reference_fai
    path model_dir, stageAs: 'models'

    output:
    tuple val(sample_id), path("${sample_id}.google_deepvariant.vcf.gz")

    script:
    def model_type = params.sample_type == "wes" ? "WES" : "WGS"
    // Extract just the basename of the model file for the relative path
    def model_basename = params.gg_deepvariant_model ? file(params.gg_deepvariant_model).name : ""
    def model_flag = model_basename ? "--customized_model models/${model_basename}" : ""
    """
    /opt/deepvariant/bin/run_deepvariant \\
        --model_type ${model_type} \\
        --ref ${reference} \\
        --reads ${bam} \\
        --output_vcf ${sample_id}.google_deepvariant.vcf.gz \\
        --num_shards ${params.threads} \\
        ${model_flag}
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
    path target_bed

    output:
    path "${sample_id}-${prefix}-*"

    script:
    def target_flag = (params.sample_type == "wes" && target_bed.name != "NO_FILE") ? "-T ${target_bed}" : ""
    """
    HGREF=${reference}
    /opt/hap.py/bin/hap.py \\
        ${benchmark_vcf} \\
        ${vcf} \\
        -f ${benchmark_bed} \\
        --threads 10 \\
        --engine vcfeval \\
        -r ${reference} \\
        ${target_flag} \\
        -o ${sample_id}-${prefix}-vcfeval
    """
}

workflow {
    // Print all parameters
    log.info ""
    log.info "========================================="
    log.info "Parabricks Variant Calling Pipeline"
    log.info "========================================="
    log.info "Sample ID         : ${params.sample_id}"
    log.info "Input             : ${params.input}"
    log.info "Reference         : ${params.reference}"
    log.info "Known Sites       : ${params.known_sites}"
    log.info "Benchmark VCF     : ${params.benchmark_vcf}"
    log.info "Benchmark BED     : ${params.benchmark_bed}"
    log.info "Target BED        : ${params.target_bed ?: 'Not provided'}"
    log.info "Sample Type       : ${params.sample_type}"
    log.info "Threads           : ${params.threads}"
    log.info "Output Directory  : ${params.output_dir}"
    log.info "========================================="
    log.info ""

    reference_ch = Channel.fromPath(params.reference).collect()
    reference_fai_ch = Channel.fromPath("${params.reference}.fai").collect()
    known_sites_ch = Channel.fromPath(params.known_sites).collect()

    // Handle optional parabricks deepvariant model file
    pb_model_ch = params.pb_deepvariant_model
        ? Channel.fromPath(params.pb_deepvariant_model).collect()
        : Channel.value(file('NO_MODEL_FILE'))

    // Handle optional google deepvariant model file
    // If model path is provided, stage the entire parent directory
    gg_model_ch = params.gg_deepvariant_model
        ? Channel.fromPath(file(params.gg_deepvariant_model).parent, type: 'dir')
        : Channel.value(file('NO_MODEL_FILE'))

    // Handle input based on input_type
    if (params.input_type == "bam") {
        // Parse BAM and BAI files
        def input_files = params.input.split(',').collect { file(it.trim()) }
        bam_ch = Channel.of([params.sample_id, input_files[0], input_files[1]])

        // Run BAM processing pipeline
        bamsort(bam_ch, reference_ch)
        markdup(bamsort.out, reference_ch)
        bqsr(markdup.out, reference_ch, known_sites_ch)

        // Extract BAM from markdup and recal file from bqsr for haplotypecaller
        markdup_bam_ch = markdup.out.map { sample_id, bam, bai -> [sample_id, bam] }
        recal_file_ch = bqsr.out.map { sample_id, bam, recal -> recal }

        haplotypecaller(markdup_bam_ch, reference_ch, recal_file_ch)

        // Run deepvariant on markdup BAM (wait for haplotypecaller to complete)
        markdup_for_dv = markdup.out.combine(haplotypecaller.out.map { it[0] })
            .map { sample_id, bam, bai, hc_sample_id -> [sample_id, bam] }
        deepvariant(markdup_for_dv, reference_ch, pb_model_ch)

        // Run google deepvariant on markdup BAM with index
        google_deepvariant(markdup.out, reference_ch, reference_fai_ch, gg_model_ch)

        // Set vcf channels for benchmarking
        haplotypecaller_vcf = haplotypecaller.out
    } else {
        // Parse FASTQ files: R1_1,R2_1,R1_2,R2_2,...
        def fq_list = params.input.split(',').collect { file(it.trim()) }
        fq_ch = Channel.of([params.sample_id, fq_list])

        germline(fq_ch, reference_ch, known_sites_ch)

        // Extract BAM only for parabricks deepvariant
        germline_bam_only = germline.out.bam.map { sample_id, bam, bai -> [sample_id, bam] }
        deepvariant(germline_bam_only, reference_ch, pb_model_ch)

        // Pass BAM with index to google deepvariant
        google_deepvariant(germline.out.bam, reference_ch, reference_fai_ch, gg_model_ch)

        // Set vcf channels for benchmarking
        haplotypecaller_vcf = germline.out.vcf.map { sample_id, vcf -> [sample_id, vcf] }
    }

    // Only run benchmark if both benchmark_vcf and benchmark_bed are provided
    if (params.benchmark_vcf && params.benchmark_bed) {
        benchmark_vcf_ch = Channel.fromPath(params.benchmark_vcf).collect()
        benchmark_vcf_tbi_ch = Channel.fromPath("${params.benchmark_vcf}.tbi").collect()
        benchmark_bed_ch = Channel.fromPath(params.benchmark_bed).collect()
        target_bed_ch = params.target_bed ? Channel.fromPath(params.target_bed).collect() : Channel.value(file("NO_FILE"))

        // Combine VCF channels with prefixes for benchmarking
        // Use haplotypecaller_vcf for germline/haplotypecaller (works for both bam and fastq input)
        vcf_prefix = params.input_type == "bam" ? "haplotypecaller" : "germline"
        vcf_ch = haplotypecaller_vcf
            .map { sample_id, vcf -> [sample_id, vcf, vcf_prefix] }
            .mix(deepvariant.out.map { sample_id, vcf -> [sample_id, vcf, "deepvariant"] })
            .mix(google_deepvariant.out.map { sample_id, vcf -> [sample_id, vcf, "google_deepvariant"] })

        benchmark(vcf_ch, reference_ch, reference_fai_ch, benchmark_vcf_ch, benchmark_vcf_tbi_ch, benchmark_bed_ch, target_bed_ch)
    }
}
