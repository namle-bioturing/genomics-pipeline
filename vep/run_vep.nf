#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Parameters
params.sample = 'PAW81754'
params.vcf = 'PAW81754.vcf.gz'
params.threads = 20
params.fasta = "/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
params.synonyms = "/references/synonyms.txt"
params.bed = "/references/hg38_Twist_Bioscience_for_Illumina_Exome_2_5_Mito.bed"
params.vep_dir = "/opt/vep/.vep/"
params.vep_cache = "/opt/vep/.vep/"
params.vep_plugins = "/opt/vep/.vep/plugins"
params.vep_fasta = "/opt/vep/.vep/homo_sapiens_merged/115_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
params.revel_plugin = "/opt/vep/.vep/plugins_data/REVEL/revel_hg38.tsv.gz"
params.alphamissense_plugin = "/opt/vep/.vep/plugins_data/AlphaMissense/AlphaMissense_hg38.tsv.gz"
params.cadd_snv_plugin = "/opt/vep/.vep/plugins_data/CADD/whole_genome_SNVs_v1.7_hg38.tsv.gz"
params.cadd_indel_plugin = "/opt/vep/.vep/plugins_data/CADD/gnomad.genomes.r4.0.indel.v1.7.hg38.tsv.gz"
params.splice_ai_snv_plugin = "/opt/vep/.vep/plugins_data/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz"
params.splice_ai_indel_plugin = "/opt/vep/.vep/plugins_data/SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz"

params.outdir = "."

// Container configuration
params.container = 'nam.le_bioturing.com/ensembl-vep:release_115.2'
params.containerOptions = "-v /mnt/nasdev2/namle/references:/references -v /mnt/nasdev2/namle/.vep:/opt/vep/.vep"

// Workflow
workflow {
    // Create input channel
    vcf_ch = Channel.fromPath(params.vcf, checkIfExists: true)

    // Run normalization
    NORMALIZE_VCF(vcf_ch, params.sample)

    // Run filtering/preprocessing
    FILTER_VCF(NORMALIZE_VCF.out.normalized_vcf, NORMALIZE_VCF.out.normalized_vcf_index, params.sample)

    // Run VEP on entire filtered VCF
    RUN_VEP(NORMALIZE_VCF.out.normalized_vcf, NORMALIZE_VCF.out.normalized_vcf_index, params.sample)
}

// Process: Normalize VCF
process NORMALIZE_VCF {
    container params.container
    containerOptions params.containerOptions
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.normalized.vcf.gz*'

    cpus params.threads
    memory '4 GB'

    input:
    path vcf
    val sample

    output:
    path "${sample}.out.normalized.vcf.gz", emit: normalized_vcf
    path "${sample}.out.normalized.vcf.gz.tbi", emit: normalized_vcf_index

    script:
    """
    echo "[NORMALIZE] Starting normalization at \$(date)"

    bcftools norm -m-any --force --check-ref w \
        -f ${params.fasta} \
        ${vcf} \
        --threads ${task.cpus} \
        -Oz -o ${sample}.out.normalized.vcf.gz

    tabix -f ${sample}.out.normalized.vcf.gz

    echo "[NORMALIZE] Normalization completed at \$(date)"
    """
}

// Process: Filter VCF with BED file
process FILTER_VCF {
    container params.container
    containerOptions params.containerOptions
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.filtered.vcf.gz*'

    cpus params.threads
    memory '4 GB'

    input:
    path normalized_vcf
    path normalized_vcf_index
    val sample

    output:
    path "${sample}.out.filtered.vcf.gz", emit: filtered_vcf
    path "${sample}.out.filtered.vcf.gz.tbi", emit: filtered_vcf_index

    script:
    """
    echo "[FILTER] Starting filtering at \$(date)"

    # Apply bed file
    bedtools intersect -header -a ${normalized_vcf} -b ${params.bed} > ${sample}.out.filtered.vcf

    # Compress and index
    bgzip -@ ${task.cpus} -f ${sample}.out.filtered.vcf
    tabix -f ${sample}.out.filtered.vcf.gz

    echo "[FILTER] Filtering completed at \$(date)"
    """
}

// Process: Run VEP on entire VCF
process RUN_VEP {
    container params.container
    containerOptions params.containerOptions
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.vep.vcf.gz*'

    cpus params.threads
    memory '20 GB'

    tag "${sample}"

    input:
    path filtered_vcf
    path filtered_vcf_index
    val sample

    output:
    path "${sample}.out.vep.vcf.gz", emit: vep_output
    path "${sample}.out.vep.vcf.gz.tbi", emit: vep_output_index

    script:
    """
    echo "[VEP] Starting VEP annotation at \$(date)"

    vep \
        --dir ${params.vep_dir} \
        --dir_cache ${params.vep_cache} \
        --dir_plugins ${params.vep_plugins} \
        --species homo_sapiens \
        --vcf \
        --offline \
        --sift b \
        --polyphen b \
        --ccds \
        --hgvs \
        --symbol \
        --numbers \
        --canonical \
        --protein \
        --af \
        --af_1kg \
        --af_gnomade \
        --af_gnomadg \
        --max_af \
        --pubmed \
        --mane \
        --tsl \
        --appris \
        --variant_class \
        --gene_phenotype \
        --mirna \
        --merged \
        --force_overwrite \
        --total_length \
        --no_stats \
        --exclude_predicted \
        --assembly GRCh38 \
        --buffer_size 20000 \
        --fork ${task.cpus} \
        --flag_pick_allele \
        --synonyms ${params.synonyms} \
        --fasta ${params.vep_fasta} \
        --plugin REVEL,file=${params.revel_plugin},no_match=1 \
        --plugin AlphaMissense,file=${params.alphamissense_plugin} \
        --plugin CADD,${params.cadd_snv_plugin},${params.cadd_indel_plugin} \
        --plugin SpliceAI,snv=${params.splice_ai_snv_plugin},indel=${params.splice_ai_indel_plugin} \
        --compress_output bgzip \
        -i ${filtered_vcf} \
        -o ${sample}.out.vep.vcf.gz

    tabix -f ${sample}.out.vep.vcf.gz

    echo "[VEP] VEP annotation completed at \$(date)"
    """
}

// nextflow run run_vep.nf --sample PAW81754 --vcf PAW81754.vcf.gz -with-timeline timeline.html -with-report report.html -with-trace trace.txt -profile standard
