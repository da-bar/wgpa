process UCSC_LAVTOPSL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-lavtopsl:482--h0b57e2e_0' :
        'quay.io/biocontainers/ucsc-lavtopsl:482--h0b57e2e_0' }"

    input:
    // One or more LAV files (e.g. the chunks of one LASTZ pair), plain or gzipped, converted in the order given.
    // They are staged in numbered folders, because the chunks of a pair have the same file name.
    tuple val(meta), path(lav, stageAs: 'lav*/*')

    output:
    tuple val(meta), path("*.psl.gz"), emit: psl
    tuple val("${task.process}"), val('ucsc'), val('482'), topic: versions, emit: versions_ucsc
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // The LAV files are joined in the order of the input list and read by lavToPsl from stdin (it reads several LAV
    // files joined into one). It writes '##' comment lines (aligner, matrix, gap penalties) for each LAV file; they
    // are removed, so the PSL holds records only.
    // The match and mismatch counts of the PSL are estimated from the percent identity of each LAV block;
    // axtChain scores the blocks from the sequences and does not use them.
    """
    set -o pipefail

    gzip -cdf ${lav instanceof List ? lav.join(' ') : lav} \\
        | lavToPsl \\
            $args \\
            stdin \\
            stdout \\
        | sed '/^#/d' \\
        | gzip -n > ${prefix}.psl.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '' | gzip -n > ${prefix}.psl.gz
    """
}
