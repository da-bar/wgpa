process LASTZ {
    tag "${meta.id} ${chunk.index + 1}/${chunk.count}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/lastz:1.04.60--hab16a5f_0' :
        'quay.io/biocontainers/lastz:1.04.60--hab16a5f_0' }"

    input:
    // chunk: [ index: 0-based chunk number, count: number of chunks of the pair, seqs: the target sequence names of this
    // chunk ]. The query sizes file lists the query sequences (twoBitInfo order). Target and query are staged in
    // separate folders, so that a self-alignment (the same 2bit file twice) does not collide.
    tuple val(meta), val(chunk), path(target_twobit, stageAs: 'target/*'), path(query_twobit, stageAs: 'query/*'), path(query_sizes, stageAs: 'query/*')
    path(matrix)

    output:
    tuple val(meta), val(chunk), path("*.lav.gz"), emit: lav
    tuple val("${task.process}"), val('lastz'), eval("lastz --version | sed -n '1s/^lastz (version \\([^ ]*\\) .*/\\1/p'"), topic: versions, emit: versions_lastz

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def matrix_arg = matrix ? "Q=${matrix}" : ''
    // Genome names are unique per FASTA file, so the same 2bit file name on both sides is a self-alignment.
    // (The name of a staged input includes its stageAs folder, e.g. 'target/x.2bit': compare the last component.)
    def self_alignment = target_twobit.toString().tokenize('/')[-1] == query_twobit.toString().tokenize('/')[-1]
    def target_seqs = chunk.seqs.join('\n')
    // As CNEr lastz(): one lastz run per (target sequence, query sequence), with the BLASTZ-style options of the
    // preset (ext.args) and Q=<matrix>; a sequence against itself runs with '--self --nomirror'. Each run writes
    // LAV ending with '# lastz end-of-file' (--markend), checked as CNEr validateLastz() does. This task runs the
    // target sequences of one chunk against every query sequence, up to task.cpus runs at a time, and joins their
    // LAV in a fixed order (target sequence, then query sequence). The LAV goes to stdout, so that no file name is
    // recorded in it: the joined LAV of all chunks is the same whatever the chunks and the number of cpus.
    """
    set -euo pipefail

    cat <<'END_TARGET_SEQS' > target_seqs.txt
${target_seqs}
END_TARGET_SEQS

    awk 'NR == FNR { q[++n] = \$1; next } { for (i = 1; i <= n; i++) print \$1 "\\t" q[i] }' \\
        $query_sizes target_seqs.txt > pairs.tsv

    run_lastz() {
        local out="lav/\$1.lav" t="\$2" q="\$3"
        if [ "${self_alignment}" = true ] && [ "\$t" = "\$q" ]; then
            lastz "${target_twobit}/\$t" --self --nomirror $args $matrix_arg --format=lav --markend > "\$out"
        else
            lastz "${target_twobit}/\$t" "${query_twobit}/\$q" $args $matrix_arg --format=lav --markend > "\$out"
        fi
    }

    mkdir lav
    n=0
    while IFS=\$'\\t' read -r t q; do
        n=\$((n + 1))
        run_lastz "\$(printf '%09d' "\$n")" "\$t" "\$q" &
        while [ "\$(jobs -rp | wc -l)" -ge ${task.cpus} ]; do wait -n; done
    done < pairs.tsv
    while [ "\$(jobs -rp | wc -l)" -gt 0 ]; do wait -n; done

    # Every run must have finished its LAV
    for i in \$(seq 1 "\$n"); do
        lav="lav/\$(printf '%09d' "\$i").lav"
        if [ "\$(tail -n 1 "\$lav" 2>/dev/null)" != '# lastz end-of-file' ]; then
            echo "ERROR: lastz run \$i (\$(sed -n "\${i}p" pairs.tsv | tr '\\t' ' ')) did not complete \$lav" >&2
            exit 1
        fi
    done

    find lav -name '*.lav' | sort | xargs cat | gzip -n > ${prefix}.lav.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '' | gzip -n > ${prefix}.lav.gz
    """
}
