# Snakefile for QIIME 2 analysis

SAMPLES = "/directory/manifest.csv"
CLASSIFIER = "/directory/gg-13-8-99-nb-classifier-qiime2019-4.qza"
METADATA = "/directory/metadata.txt"
OUTPUT_DIR = "/directory"
NUM_THREADS = 8

rule all:
    input:
        f"{OUTPUT_DIR}/feature-table-taxonomy.biom",
        f"{OUTPUT_DIR}/alpha_rarefaction.qzv",
        f"{OUTPUT_DIR}/core-metrics-results/faith_pd_vector.qza",
        f"{OUTPUT_DIR}/core-metrics-results/shannon_vector.qza",
        f"{OUTPUT_DIR}/taxonomy.txt"

rule import_data:
    output:
        f"{OUTPUT_DIR}/paired-end-demux.qza"
    shell:
        "qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' "
        "--input-path {SAMPLES} --output-path {output} "
        "--input-format PairedEndFastqManifestPhred33"

rule summarize_demux:
    input:
        f"{OUTPUT_DIR}/paired-end-demux.qza"
    output:
        f"{OUTPUT_DIR}/paired-end-demux.qzv"
    shell:
        "qiime demux summarize --i-data {input} --o-visualization {output}"

rule dada2_denoise:
    input:
        f"{OUTPUT_DIR}/paired-end-demux.qza"
    output:
        table=f"{OUTPUT_DIR}/table.qza",
        rep_seqs=f"{OUTPUT_DIR}/rep-seqs.qza",
        stats=f"{OUTPUT_DIR}/dada2denoising-stats.qza"
    shell:
        "qiime dada2 denoise-paired --i-demultiplexed-seqs {input} "
        "--p-trunc-len-f 270 --p-trunc-len-r 270 "
        "--p-trim-left-f 20 --p-trim-left-r 20 "
        "--p-max-ee 6 --p-n-threads {NUM_THREADS} "
        "--o-table {output.table} --o-representative-sequences {output.rep_seqs} "
        "--o-denoising-stats {output.stats}"

rule phylogeny:
    input:
        f"{OUTPUT_DIR}/rep-seqs.qza"
    output:
        alignment=f"{OUTPUT_DIR}/aligned_rep_seqs_dada2.qza",
        masked=f"{OUTPUT_DIR}/masked_aligned_rep_seqs_dada2.qza",
        unrooted=f"{OUTPUT_DIR}/unrooted_tree.qza",
        rooted=f"{OUTPUT_DIR}/rooted_tree.qza"
    shell:
        "qiime phylogeny align-to-tree-mafft-fasttree --i-sequences {input} "
        "--o-alignment {output.alignment} --o-masked-alignment {output.masked} "
        "--o-tree {output.unrooted} --o-rooted-tree {output.rooted}"

rule classify_taxonomy:
    input:
        rep_seqs=f"{OUTPUT_DIR}/rep-seqs.qza"
    output:
        f"{OUTPUT_DIR}/taxonomy.qza"
    shell:
        "qiime feature-classifier classify-sklearn --i-classifier {CLASSIFIER} "
        "--i-reads {input.rep_seqs} --o-classification {output}"

rule filter_samples:
    input:
        table=f"{OUTPUT_DIR}/table.qza",
        taxonomy=f"{OUTPUT_DIR}/taxonomy.qza"
    output:
        f"{OUTPUT_DIR}/filtered-table.qza"
    shell:
        "qiime taxa filter-table --i-table {input.table} "
        "--i-taxonomy {input.taxonomy} --p-exclude mitochondria,chloroplast "
        "--o-filtered-table {output}"

rule alpha_rarefaction:
    input:
        table=f"{OUTPUT_DIR}/filtered-table.qza",
        phylogeny=f"{OUTPUT_DIR}/rooted_tree.qza"
    output:
        f"{OUTPUT_DIR}/alpha_rarefaction.qzv"
    shell:
        "qiime diversity alpha-rarefaction --i-table {input.table} "
        "--i-phylogeny {input.phylogeny} --o-visualization {output} "
        "--p-max-depth 18000"

rule core_metrics:
    input:
        phylogeny=f"{OUTPUT_DIR}/rooted_tree.qza",
        table=f"{OUTPUT_DIR}/filtered-table.qza"
    output:
        directory(f"{OUTPUT_DIR}/core-metrics-results")
    shell:
        "qiime diversity core-metrics-phylogenetic --i-phylogeny {input.phylogeny} "
        "--i-table {input.table} --p-sampling-depth 10000 "
        "--m-metadata-file {METADATA} --output-dir {output}"

rule extract_distance_matrices:
    input:
        weighted=f"{OUTPUT_DIR}/core-metrics-results/weighted_unifrac_distance_matrix.qza",
        unweighted=f"{OUTPUT_DIR}/core-metrics-results/unweighted_unifrac_distance_matrix.qza",
        bray_curtis=f"{OUTPUT_DIR}/core-metrics-results/bray_curtis_distance_matrix.qza"
    output:
        directory(f"{OUTPUT_DIR}/core-metrics-results/")
    shell:
        "qiime tools extract --input-path {input.weighted} --output-path {output} && "
        "qiime tools extract --input-path {input.unweighted} --output-path {output} && "
        "qiime tools extract --input-path {input.bray_curtis} --output-path {output}"

rule export_alpha_metrics:
    input:
        faith=f"{OUTPUT_DIR}/core-metrics-results/faith_pd_vector.qza",
        shannon=f"{OUTPUT_DIR}/core-metrics-results/shannon_vector.qza"
    output:
        directory(f"{OUTPUT_DIR}/core-metrics-results/")
    shell:
        "qiime tools extract --input-path {input.faith} --output-path {output} && "
        "qiime tools extract --input-path {input.shannon} --output-path {output}"

rule export_taxonomy:
    input:
        f"{OUTPUT_DIR}/taxonomy.qza"
    output:
        f"{OUTPUT_DIR}/taxonomy.txt"
    shell:
        "qiime tools export --input-path {input} --output-path {output}"

rule export_feature_table:
    input:
        f"{OUTPUT_DIR}/core-metrics-results/rarefied_table.qza"
    output:
        f"{OUTPUT_DIR}/rarefied_table.biom"
    shell:
        "qiime tools extract --input-path {input} --output-path {output}"

rule merge_taxonomy_with_table:
    input:
        biom=f"{OUTPUT_DIR}/feature-table.biom",
        taxonomy=f"{OUTPUT_DIR}/taxonomy.txt"
    output:
        f"{OUTPUT_DIR}/feature-table-taxonomy.biom"
    shell:
        "biom add-metadata -i {input.biom} -o {output} "
        "--observation-metadata-fp {input.taxonomy} --sc-separated taxonomy"
