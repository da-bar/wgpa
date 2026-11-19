//
// Align each query genome to its target genome and convert the alignments to PSL, with the pair's aligner:
// LAST (the default), LASTZ or MUMmer4 (meta.aligner 'lastz' or 'mummer'). A pair with an alignment in the
// samplesheet (column alignment: MAF or PSL, optionally gzipped) is not aligned: its alignment is checked
// against the pair's genomes, and a MAF is converted to PSL.
//
//   LAST:    lastdb -c, lastal (CNEr options and matrix), maf-convert psl
//   LASTZ:   lastz for each (target sequence, query sequence) with the CNEr options and matrix, lavToPsl. One task per
//            chunk of target sequences (a sequence of at least lastz_chunk_size bases, or consecutive smaller ones up
//            to that size together); the chunks are joined back in target order, so the PSL does not depend on them.
//   MUMmer4: nucmer, then delta to PSL (CONVERT_DELTA2PSL)
//

include { LAST_LASTDB       } from '../../../modules/nf-core/last/lastdb'
include { LAST_LASTAL       } from '../../../modules/nf-core/last/lastal'
include { LAST_MAFCONVERT   } from '../../../modules/nf-core/last/mafconvert'
include { LASTZ             } from '../../../modules/local/lastz'
include { UCSC_LAVTOPSL     } from '../../../modules/local/ucsc/lavtopsl'
include { MUMMER4_NUCMER    } from '../../../modules/local/mummer4/nucmer'
include { CONVERT_DELTA2PSL } from '../../../modules/local/convert/delta2psl'
include { CHECK_ALIGNMENT   } from '../../../modules/local/check/alignment'

workflow PAIRWISE_ALIGN {

    take:
    ch_samplesheet   // channel: [ val(meta), path(target_fasta), path(query_fasta), path(alignment) ], alignment = [] if not given; meta.alignment_format = 'maf' or 'psl' if given
    ch_lastal_matrix // channel: [ val(meta), path(lastal_matrix) ], the preset's scoring matrix for lastal -p
    ch_lastz_matrix  // channel: [ val(meta), path(lastz_matrix) ], the preset's scoring matrix for lastz Q= (CNEr: the axtChain matrix)
    ch_genomes       // channel: [ val(meta), path(target_2bit), path(target_sizes), path(query_2bit), path(query_sizes) ], each pair's genomes
    lastz_chunk_size // integer: target bases per LASTZ task (a longer target sequence is a chunk of its own)

    main:

    def ch_sizes = ch_genomes
        .map { meta, _target_twobit, target_sizes, _query_twobit, query_sizes -> [ meta, target_sizes, query_sizes ] }

    //
    // Pairs with an alignment input are not aligned. The other pairs go to their aligner (meta.aligner is only
    // set for LASTZ and MUMmer4, so that the meta map of a LAST pair is as before).
    //
    def ch_input = ch_samplesheet
        .branch { meta, _target, _query, _alignment ->
            alignment: meta.alignment_format != null
            lastz: meta.aligner == 'lastz'
            mummer: meta.aligner == 'mummer'
            lastal: true
        }
    def ch_pairs = ch_input.lastal
        .map { meta, target, query, _alignment -> [ meta, target, query ] }

    //
    // Build one LAST index per distinct (target FASTA, lastdb arguments)
    //
    def ch_lastdb_in = ch_pairs
        .map { meta, target, _query ->
            [ [ id: meta.target_name, lastdb_args: lastdbArgs(meta.preset) ], target ]
        }
        .unique { meta, target -> [ target.toUriString(), meta.lastdb_args ] }

    LAST_LASTDB ( ch_lastdb_in )

    //
    // LAST_LASTAL takes the query and the index as two separate inputs. Build one combined channel per
    // pair and split it with multiMap so that the two inputs stay in lockstep.
    // The target name is unique per FASTA file (checked in PIPELINE_INITIALISATION).
    //
    def ch_index = LAST_LASTDB.out.index
        .map { meta, index -> [ [ meta.id, meta.lastdb_args ], index ] }

    def ch_lastal_in = ch_pairs
        .join(ch_lastal_matrix, failOnDuplicate: true)
        .map { meta, _target, query, matrix -> [ [ meta.target_name, lastdbArgs(meta.preset) ], meta, query, matrix ] }
        .combine(ch_index, by: 0)
        .multiMap { _key, meta, query, matrix, index ->
            query: [ meta, query, matrix ]
            index: index
        }

    LAST_LASTAL ( ch_lastal_in.query, ch_lastal_in.index )

    //
    // LASTZ on the 2bit files (as CNEr): one task per chunk of target sequences, each running one lastz per
    // (target sequence, query sequence). The LAV of the chunks is joined in chunk order (= target order) and
    // converted to PSL, so the PSL is the same whatever the chunk size.
    //
    def ch_lastz_in = ch_input.lastz
        .map { meta, _target, _query, _alignment -> [ meta ] }
        .join(ch_genomes, failOnDuplicate: true)
        .join(ch_lastz_matrix, failOnDuplicate: true)
        .flatMap { meta, target_twobit, target_sizes, query_twobit, query_sizes, matrix ->
            def chunks = sequenceChunks(target_sizes, lastz_chunk_size)
            chunks.withIndex().collect { seqs, index ->
                [ meta, [ index: index, count: chunks.size(), seqs: seqs ], target_twobit, query_twobit, query_sizes, matrix ]
            }
        }
        .multiMap { meta, chunk, target_twobit, query_twobit, query_sizes, matrix ->
            genomes: [ meta, chunk, target_twobit, query_twobit, query_sizes ]
            matrix: matrix
        }

    LASTZ ( ch_lastz_in.genomes, ch_lastz_in.matrix )

    // Gather the chunks of each pair as soon as all of them are done, in chunk order
    def ch_lastz_lav = LASTZ.out.lav
        .map { meta, chunk, lav -> [ groupKey(meta, chunk.count), [ chunk.index, lav ] ] }
        .groupTuple()
        .map { key, chunks -> [ key.getGroupTarget(), chunks.sort { chunk -> chunk[0] }.collect { chunk -> chunk[1] } ] }

    UCSC_LAVTOPSL ( ch_lastz_lav )

    //
    // MUMmer4 nucmer on the FASTA files, then delta to PSL
    //
    MUMMER4_NUCMER (
        ch_input.mummer.map { meta, target, query, _alignment -> [ meta, target, query ] }
    )
    CONVERT_DELTA2PSL ( MUMMER4_NUCMER.out.delta )

    //
    // Alignment inputs: check that their target and query sequences are those of the pair's genomes
    // (names and lengths), so that a wrong genome or swapped target and query stop the pipeline here
    //
    CHECK_ALIGNMENT (
        ch_input.alignment
            .map { meta, _target, _query, alignment -> [ meta, alignment ] }
            .join(ch_sizes, failOnDuplicate: true)
    )
    def ch_checked = CHECK_ALIGNMENT.out.alignment
        .branch { meta, _alignment ->
            maf: meta.alignment_format == 'maf'
            psl: true
        }

    //
    // MAF to PSL, the input of axtChain (no reference files are needed for PSL). A PSL input goes to the
    // chaining as it is.
    //
    LAST_MAFCONVERT (
        LAST_LASTAL.out.maf
            .mix(ch_checked.maf)
            .map { meta, maf -> [ meta, maf, 'psl' ] },
        [ [:], [], [], [], [], [] ]
    )

    emit:
    maf   = LAST_LASTAL.out.maf                     // channel: [ val(meta), path(maf.gz) ], pairs aligned with LAST
    stats = LAST_LASTAL.out.multiqc                 // channel: [ val(meta), path(tsv) ], pairs aligned with LAST
    delta = MUMMER4_NUCMER.out.delta                // channel: [ val(meta), path(delta.gz) ], pairs aligned with MUMmer4
    psl   = LAST_MAFCONVERT.out.alignment           // channel: [ val(meta), path(psl) ], every pair
        .mix(ch_checked.psl, UCSC_LAVTOPSL.out.psl, CONVERT_DELTA2PSL.out.psl)
}

//
// lastdb arguments for a preset. CNEr indexes the target with 'lastdb -c' for every preset (lowercase,
// i.e. soft-masked, bases are excluded from initial matches), so pairs with the same target FASTA share
// one index whatever their preset. The value is passed to LAST_LASTDB as meta.lastdb_args
// (conf/modules.config) because it is part of the key that decides which pairs share an index.
//
def lastdbArgs(preset) {
    def args = [
        near   : '-c',
        medium : '-c',
        far    : '-c',
    ]
    if (!args.containsKey(preset)) {
        error("Unknown preset '${preset}': expected one of ${args.keySet().join(', ')}")
    }
    return args[preset]
}

//
// Split the sequences of a genome (its sizes file, in 2bit order) into chunks, for the LASTZ tasks (target) and the
// chaining tiles (CHAIN_NET, target and query): a sequence of at least chunk_size bases is a chunk of its own;
// consecutive shorter sequences go together while their total stays within chunk_size. The chunks keep the order of
// the sequences.
//
def sequenceChunks(sizes, chunk_size) {
    def chunks = []
    def current = []
    def current_size = 0L
    sizes.readLines().findAll { line -> line.trim() }.each { line ->
        def fields = line.split('\t')
        def name = fields[0]
        def size = fields[1] as long
        if (current && current_size + size > chunk_size) {
            chunks << current
            current = []
            current_size = 0L
        }
        current << name
        current_size += size
    }
    if (current) {
        chunks << current
    }
    return chunks
}
