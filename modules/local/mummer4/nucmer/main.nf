process MUMMER4_NUCMER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mummer4:4.0.1--pl5321h9948957_0' :
        'quay.io/biocontainers/mummer4:4.0.1--pl5321h9948957_0' }"

    input:
    // Target (reference) and query are staged in separate folders, so that a self-alignment (the same FASTA twice) does not collide
    tuple val(meta), path(target_fasta, stageAs: 'target/*'), path(query_fasta, stageAs: 'query/*')

    output:
    tuple val(meta), path("*.delta.gz"), emit: delta
    tuple val("${task.process}"), val('mummer4'), eval("nucmer --version"), topic: versions, emit: versions_mummer4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // nucmer runs with task.cpus threads. With more than one thread it finds the same alignments, but writes them
    // in a run-dependent order; CONVERT_DELTA2PSL sorts them, so the PSL and everything made from it do not depend
    // on the number of threads (the delta file does). Gzipped FASTA is decompressed first. nucmer writes the
    // absolute paths of the two FASTA files in the first line of the delta file; they are replaced by the
    // input file names, so the delta file does not depend on the work directory.
    // The name of a staged input includes its stageAs folder ('target/x.fa'): use the last component
    def target_file  = target_fasta.toString().tokenize('/')[-1]
    def query_file   = query_fasta.toString().tokenize('/')[-1]
    def target_plain = target_file.endsWith('.gz') ? "target/plain/${target_file - ~/\.gz$/}" : "${target_fasta}"
    def query_plain  = query_file.endsWith('.gz')  ? "query/plain/${query_file - ~/\.gz$/}"   : "${query_fasta}"
    """
    set -euo pipefail

    if [ "${target_plain}" != "${target_fasta}" ]; then
        mkdir -p target/plain
        gzip -cd ${target_fasta} > ${target_plain}
    fi
    if [ "${query_plain}" != "${query_fasta}" ]; then
        mkdir -p query/plain
        gzip -cd ${query_fasta} > ${query_plain}
    fi

    nucmer \\
        --threads=${task.cpus} \\
        $args \\
        --delta=${prefix}.delta \\
        ${target_plain} \\
        ${query_plain}

    sed '1s|.*|${target_file} ${query_file}|' ${prefix}.delta | gzip -n > ${prefix}.delta.gz
    rm ${prefix}.delta
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '' | gzip -n > ${prefix}.delta.gz
    """
}
