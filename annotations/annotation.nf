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

// InterVar local paths
params.intervar_dir = "/mnt/nasdev2/namle/references/intervar"
params.annovar_dir = "/mnt/nasdev2/namle/references/annovar"
// params.intervar_py = "/mnt/nasdev2/namle/references/InterVar/Intervar.py"
// params.intervar_config = "/mnt/nasdev2/namle/references/config.ini"
// params.intervar_db = "/mnt/nasdev2/namle/references/intervardb"

// Annotate container configuration
params.annotate_container = 'nam.le_bioturing.com/annotate:latest'
params.annotate_containerOptions = "-v /mnt/nasdev2/namle/references:/references"

// Reference files for annotation
params.hgnc_mapping = "/references/hgnc.txt"
params.inheritance_genes = "/references/inheritance_genes.tsv"
params.omim_info = "/references/omim_info.txt"

// Workflow
workflow {
    // Create input channel
    vcf_ch = Channel.fromPath(params.vcf, checkIfExists: true)
    intervar_ch = Channel.fromPath(params.intervar_dir, checkIfExists: true)
    annovar_ch = Channel.fromPath(params.annovar_dir, checkIfExists: true)

    // Run normalization
    NORMALIZE_VCF(vcf_ch, params.sample)

    // Run filtering/preprocessing
    FILTER_VCF(NORMALIZE_VCF.out.normalized_vcf, NORMALIZE_VCF.out.normalized_vcf_index, params.sample)

    // Run VEP on entire filtered VCF
    RUN_VEP(NORMALIZE_VCF.out.normalized_vcf, NORMALIZE_VCF.out.normalized_vcf_index, params.sample)

    // Run InterVar on normalized VCF
    RUN_INTERVAR(NORMALIZE_VCF.out.normalized_vcf, NORMALIZE_VCF.out.normalized_vcf_index, params.sample, intervar_ch, annovar_ch)

    // Run annotation script combining VEP and InterVar results
    ANNOTATE_VCF(
        RUN_VEP.out.vep_output,
        RUN_VEP.out.vep_output_index,
        RUN_INTERVAR.out.intervar_classification,
        params.sample
    )
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
    path "${sample}.out.vep.bcf", emit: vep_output
    path "${sample}.out.vep.bcf.csi", emit: vep_output_index

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
        --flag_pick_allele_gene \
        --synonyms ${params.synonyms} \
        --fasta ${params.vep_fasta} \
        --plugin REVEL,file=${params.revel_plugin},no_match=1 \
        --plugin AlphaMissense,file=${params.alphamissense_plugin} \
        --plugin CADD,${params.cadd_snv_plugin},${params.cadd_indel_plugin} \
        --plugin SpliceAI,snv=${params.splice_ai_snv_plugin},indel=${params.splice_ai_indel_plugin} \
        -i ${filtered_vcf} \
        -o STDOUT | bcftools view -O b -o ${sample}.out.vep.bcf -

    bcftools index -f ${sample}.out.vep.bcf

    echo "[VEP] VEP annotation completed at \$(date)"
    """
}

// Process: Run InterVar annotation (runs locally, not in Docker)
process RUN_INTERVAR {
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.intervar*'

    cpus 10
    memory '20 GB'

    tag "${sample}"

    input:
    path normalized_vcf
    path normalized_vcf_index
    val sample
    path intervar_dir
    path annovar_dir

    output:
    path "${sample}.intervar.hg38_multianno.txt", emit: intervar_output
    path "${sample}.intervar.hg38_multianno.txt.intervar", emit: intervar_classification

    script:

    def intervar_script = "${intervar_dir}/Intervar.py"
    def intervar_config = "${intervar_dir}/config.ini"
    def intervar_db = "${intervar_dir}/intervardb"

    """
    echo "[INTERVAR] Starting InterVar annotation at \$(date)"
    echo "[INTERVAR] Working directory: \$(pwd)"

    # Run InterVar locally
    python ${intervar_script} \
        -c ${intervar_config} \
        -b hg38 \
        -i ${normalized_vcf} \
        --input_type=VCF \
        -o ${sample}.intervar \
        --threads ${task.cpus} \
        -t ${intervar_db}

    echo "[INTERVAR] InterVar annotation completed at \$(date)"
    """
}

// Process: Annotate VCF with combined annotations
process ANNOTATE_VCF {
    container params.annotate_container
    containerOptions params.annotate_containerOptions
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.parquet'

    cpus 5
    memory '32 GB'

    tag "${sample}"

    input:
    path vep_vcf
    path vep_vcf_index
    path intervar_file
    val sample

    output:
    path "${sample}.annotated.parquet", emit: annotated_parquet

    script:
    """
    echo "[ANNOTATE] Starting combined annotation at \$(date)"

    # Run annotate.py with all required inputs
    annotate.py \
        --vcf ${vep_vcf} \
        --hgnc-mapping ${params.hgnc_mapping} \
        --inheritance-genes ${params.inheritance_genes} \
        --omim-info ${params.omim_info} \
        --intervar-file ${intervar_file} \
        --output ${sample}.annotated.parquet \
        --processes ${task.cpus}

    echo "[ANNOTATE] Combined annotation completed at \$(date)"
    echo "[ANNOTATE] Output file: ${sample}.annotated.parquet"
    """
}

// nextflow run annotation.nf -resume --sample PAW81754 --vcf PAW81754.vcf.gz -with-timeline timeline.html -with-report report.html -with-trace trace.txt -profile standard
