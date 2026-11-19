# Convert a MUMmer delta file (NUCMER) to PSL, one record per delta alignment.
#
# Delta alignment header: rs re qs qe errors simerrors stops (1-based, inclusive; qs > qe on the reverse strand,
# where the query coordinates are those of the forward strand). Then one integer per indel, ending with 0:
# d > 0: d-1 aligned pairs, then one target (reference) base aligned to a gap in the query;
# d < 0: -d-1 aligned pairs, then one query base aligned to a gap in the target.
# errors counts mismatches plus gap bases, so mismatches = errors - gap bases.
function flush_block() {
    if (blen > 0) {
        n_blocks++
        sizes = sizes blen ","
        q_starts = q_starts bq ","
        t_starts = t_starts bt ","
        aligned += blen
    }
}
function new_block(tpos, qpos) { blen = 0; bt = tpos; bq = qpos }
BEGIN { OFS = "\t" }
NR == 1 { next }
NR == 2 {
    if ($1 != "NUCMER") { print "ERROR: not a NUCMER delta file (line 2 is '" $0 "')" > "/dev/stderr"; failed = 1; exit 1 }
    next
}
/^>/ {
    t_name = substr($1, 2); q_name = $2; t_size = $3; q_size = $4
    in_aln = 0
    next
}
!in_aln && NF == 7 {
    rs = $1; re = $2; qs = $3; qe = $4; errors = $5
    strand = qs <= qe ? "+" : "-"
    # 0-based start of the alignment on the target and on the query strand that is aligned
    t = rs - 1
    q = strand == "+" ? qs - 1 : q_size - qs
    t_end = re; q_end = strand == "+" ? qe : q_size - qe + 1
    n_blocks = 0; sizes = ""; q_starts = ""; t_starts = ""; aligned = 0
    t_ins_n = 0; t_ins_b = 0; q_ins_n = 0; q_ins_b = 0; last_gap = ""
    new_block(t, q)
    in_aln = 1
    next
}
in_aln && NF == 1 {
    d = $1 + 0
    if (d == 0) {
        # The rest of the alignment is gapless
        blen += t_end - t
        q += t_end - t
        t = t_end
        flush_block()
        if (q != q_end) { print "ERROR: alignment " t_name " " rs "-" re " / " q_name " " qs "-" qe " does not add up (query " q " != " q_end ")" > "/dev/stderr"; failed = 1; exit 1 }
        gap_bases = t_ins_b + q_ins_b
        mismatches = errors - gap_bases
        if (mismatches < 0) { print "ERROR: alignment " t_name " " rs "-" re " / " q_name " " qs "-" qe " has fewer errors than gap bases" > "/dev/stderr"; failed = 1; exit 1 }
        q_start = strand == "+" ? qs - 1 : qe - 1
        q_stop  = strand == "+" ? qe : qs
        print aligned - mismatches, mismatches, 0, 0, q_ins_n, q_ins_b, t_ins_n, t_ins_b, strand, \
              q_name, q_size, q_start, q_stop, t_name, t_size, rs - 1, re, n_blocks, sizes, q_starts, t_starts
        in_aln = 0
        next
    }
    step = (d > 0 ? d : -d) - 1
    blen += step; t += step; q += step
    flush_block()
    # Consecutive gap bases of one type (d = 1 or -1 after a gap) are one insert; an aligned pair ends it
    if (step > 0) last_gap = ""
    if (d > 0) {
        # target base against a gap in the query
        if (last_gap != "t") t_ins_n++
        t_ins_b++; t++; last_gap = "t"
    } else {
        if (last_gap != "q") q_ins_n++
        q_ins_b++; q++; last_gap = "q"
    }
    new_block(t, q)
    next
}
{ print "ERROR: line " NR " of the delta file is not valid: '" $0 "'" > "/dev/stderr"; failed = 1; exit 1 }
END { if (!failed && in_aln) { print "ERROR: the delta file ends inside an alignment" > "/dev/stderr"; exit 1 } }
