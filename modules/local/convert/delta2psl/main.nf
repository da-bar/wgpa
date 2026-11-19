process CONVERT_DELTA2PSL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a1/a125c778baf3865331101a104b60d249ee15fe1dca13bdafd888926cc5490a34/data' :
        'community.wave.seqera.io/library/gawk:5.3.1--e09efb5dfc4b8156' }"

    input:
    tuple val(meta), path(delta)

    output:
    tuple val(meta), path("*.psl.gz"), emit: psl
    tuple val("${task.process}"), val('gawk'), eval("awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//'"), topic: versions, emit: versions_gawk

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // The converter is delta2psl.awk next to this file; its text is part of the task script (and hash).
    // One PSL record per delta alignment, with exact match and mismatch counts (no repMatches or nCount).
    // The records are sorted by target, target position, query, query position and strand (then by the whole
    // line), so the PSL does not depend on the order in which nucmer wrote the alignments.
    def awk_program = file("${moduleDir}/delta2psl.awk", checkIfExists: true).text
    """
    set -euo pipefail

    cat <<'END_DELTA2PSL' > delta2psl.awk
${awk_program}
END_DELTA2PSL

    gzip -cdf $delta \\
        | awk -f delta2psl.awk \\
        | LC_ALL=C sort -t \$'\\t' -k14,14 -k16,16n -k17,17n -k10,10 -k12,12n -k13,13n -k9,9 \\
        | gzip -n > ${prefix}.psl.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '' | gzip -n > ${prefix}.psl.gz
    """
}
