process CONVERT_SPLITPSL {
    tag "$meta.id"
    label 'process_single'
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a1/a125c778baf3865331101a104b60d249ee15fe1dca13bdafd888926cc5490a34/data' :
        'community.wave.seqera.io/library/gawk:5.3.1--e09efb5dfc4b8156' }"

    input:
    // target_chunks, query_chunks: lists of chunks, each a list of sequence names
    tuple val(meta), path(psl), val(target_chunks), val(query_chunks)

    output:
    tuple val(meta), path("tiles/*.psl"), emit: psl
    tuple val("${task.process}"), val('gawk'), eval("awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//'"), topic: versions, emit: versions_gawk

    when:
    task.ext.when == null || task.ext.when

    script:
    // Splits the PSL into tiles: one file per (target chunk, query chunk) with records, tiles/t<i>_q<j>.psl (chunk
    // numbers from 1, zero-padded, so that sorting the file names sorts the tiles). A record of a sequence that is in
    // no chunk is an error. Lines that are not PSL records (a psLayout header) are left out. Without any record the
    // output is one empty tile, tiles/empty.psl. The tiles are plain text, written in one pass over the PSL.
    def width = [3, "${[target_chunks.size(), query_chunks.size()].max()}".size()].max()
    def table = { chunks -> chunks.withIndex().collectMany { seqs, i -> seqs.collect { name -> "${name}\t${String.format("%0${width}d", i + 1)}" } }.join('\n') }
    """
    set -euo pipefail

    cat <<'END_TARGET_CHUNKS' > target_chunks.tsv
${table(target_chunks)}
END_TARGET_CHUNKS
    cat <<'END_QUERY_CHUNKS' > query_chunks.tsv
${table(query_chunks)}
END_QUERY_CHUNKS

    mkdir tiles
    gzip -cdf $psl \\
        | awk -F '\\t' '
            FILENAME == "target_chunks.tsv" { t[\$1] = \$2; next }
            FILENAME == "query_chunks.tsv"  { q[\$1] = \$2; next }
            \$1 !~ /^[0-9]+\$/ { next }
            !(\$14 in t) { print "ERROR: target sequence " \$14 " of the PSL is in no chunk" > "/dev/stderr"; exit 1 }
            !(\$10 in q) { print "ERROR: query sequence " \$10 " of the PSL is in no chunk" > "/dev/stderr"; exit 1 }
            { print > ("tiles/t" t[\$14] "_q" q[\$10] ".psl") }
        ' target_chunks.tsv query_chunks.tsv -

    if [ -z "\$(ls tiles)" ]; then
        touch tiles/empty.psl
    fi
    """

    stub:
    """
    mkdir tiles
    touch tiles/empty.psl
    """
}
