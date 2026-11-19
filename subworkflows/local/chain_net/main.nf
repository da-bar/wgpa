//
// Chain and net the pairwise alignments with the UCSC kent tools, as CNEr (and v1.0.0) does
//
// With chain_chunk_size > 0 the alignments are chained in tiles, as UCSC's own pipeline does for large genomes: the
// PSL is split by (target chunk, query chunk) (sequenceChunks, CONVERT_SPLITPSL), axtChain runs on each tile and
// chainMergeSort joins the tiles. A chain links blocks of one target and one query sequence, so the tiles make the
// same chains as one axtChain over the whole PSL, with a fraction of its memory. Only the order of chains with the
// same score differs (axtChain orders them by its own internal order, chainMergeSort by tile), and so do the chain
// IDs, and the nets where equal-scoring chains compete. chain_chunk_size = 0 (the default) runs one axtChain per
// pair, as v1.0.0.
//

include { CONVERT_SPLITPSL    } from '../../../modules/local/convert/splitpsl'
include { UCSC_AXTCHAIN       } from '../../../modules/local/ucsc/axtchain'
include { UCSC_CHAINMERGESORT } from '../../../modules/local/ucsc/chainmergesort'
include { UCSC_CHAINPRENET    } from '../../../modules/local/ucsc/chainprenet'
include { UCSC_CHAINNET       } from '../../../modules/local/ucsc/chainnet'
include { UCSC_NETSYNTENIC    } from '../../../modules/local/ucsc/netsyntenic'
include { UCSC_NETTOAXT       } from '../../../modules/local/ucsc/nettoaxt'
include { UCSC_AXTSORT        } from '../../../modules/local/ucsc/axtsort'
include { sequenceChunks      } from '../pairwise_align'

workflow CHAIN_NET {

    take:
    ch_psl          // channel: [ val(meta), path(psl) ]
    ch_genomes      // channel: [ val(meta), path(target_2bit), path(target_sizes), path(query_2bit), path(query_sizes) ]
    ch_score_scheme // channel: [ val(meta), path(axtchain_matrix) ], the preset's matrix for axtChain -scoreScheme
    chain_chunk_size //   long: target and query bases per chaining tile; 0 = one axtChain per pair

    main:

    def ch_twobits = ch_genomes
        .map { meta, target_twobit, _target_sizes, query_twobit, _query_sizes -> [ meta, target_twobit, query_twobit ] }
    def ch_sizes = ch_genomes
        .map { meta, _target_twobit, target_sizes, _query_twobit, query_sizes -> [ meta, target_sizes, query_sizes ] }

    //
    // PSL -> chains, for the whole PSL or for each tile. UCSC_AXTCHAIN takes the score scheme as a separate input:
    // split one combined channel with multiMap so that each pair gets its own preset's matrix.
    //
    def ch_chain_psl
    if (chain_chunk_size > 0) {
        CONVERT_SPLITPSL (
            ch_psl
                .join(ch_sizes, failOnDuplicate: true)
                .map { meta, psl, target_sizes, query_sizes ->
                    [ meta, psl, sequenceChunks(target_sizes, chain_chunk_size), sequenceChunks(query_sizes, chain_chunk_size) ]
                }
        )
        // one item per tile; meta.chain_tile names the tile (the output file stem), meta.chain_tiles counts them
        ch_chain_psl = CONVERT_SPLITPSL.out.psl
            .flatMap { meta, tiles ->
                def tile_list = [tiles].flatten()
                tile_list.collect { tile -> [ meta + [ chain_tile: tile.baseName, chain_tiles: tile_list.size() ], tile ] }
            }
    } else {
        ch_chain_psl = ch_psl
    }

    def ch_axtchain_in = ch_chain_psl
        .map { meta, psl -> [ pairMeta(meta), meta, psl ] }
        .combine(ch_twobits.join(ch_score_scheme, failOnDuplicate: true), by: 0)
        .multiMap { _pair_meta, meta, psl, target_twobit, query_twobit, score_scheme ->
            psl: [ meta, psl, target_twobit, query_twobit ]
            score_scheme: score_scheme
        }

    UCSC_AXTCHAIN ( ch_axtchain_in.psl, ch_axtchain_in.score_scheme )

    // The chains of a pair: one file, or one per tile (gathered once all the tiles of the pair are chained)
    def ch_pair_chains = UCSC_AXTCHAIN.out.chain
        .branch { meta, _chain ->
            tiled: meta.chain_tile != null
            whole: true
        }

    // <target>_<query>.all.chain
    UCSC_CHAINMERGESORT (
        ch_pair_chains.whole.mix(
            ch_pair_chains.tiled
                .map { meta, chain -> [ groupKey(pairMeta(meta), meta.chain_tiles), chain ] }
                .groupTuple()
                .map { key, chains -> [ key.getGroupTarget(), chains ] }
        )
    )

    // <target>_<query>.all.pre.chain
    UCSC_CHAINPRENET (
        UCSC_CHAINMERGESORT.out.chain.join(ch_sizes, failOnDuplicate: true)
    )

    // <target>_<query>.target.net and .query.net
    UCSC_CHAINNET (
        UCSC_CHAINPRENET.out.chain.join(ch_sizes, failOnDuplicate: true)
    )

    // <target>_<query>.noClass.net: the target-referenced net with synteny information
    UCSC_NETSYNTENIC ( UCSC_CHAINNET.out.target_net )

    // <target>_<query>.net.axt: the alignments of the net, sorted by target position
    UCSC_NETTOAXT (
        UCSC_NETSYNTENIC.out.net
            .join(UCSC_CHAINPRENET.out.chain, failOnDuplicate: true)
            .join(ch_twobits, failOnDuplicate: true)
    )
    UCSC_AXTSORT ( UCSC_NETTOAXT.out.axt )

    emit:
    chain      = UCSC_AXTCHAIN.out.chain       // channel: [ val(meta), path(chain.gz) ], one per tile (meta.chain_tile) with chain_chunk_size > 0
    all_chain  = UCSC_CHAINMERGESORT.out.chain // channel: [ val(meta), path(all.chain.gz) ]
    pre_chain  = UCSC_CHAINPRENET.out.chain    // channel: [ val(meta), path(all.pre.chain.gz) ]
    target_net = UCSC_CHAINNET.out.target_net  // channel: [ val(meta), path(target.net.gz) ]
    query_net  = UCSC_CHAINNET.out.query_net   // channel: [ val(meta), path(query.net.gz) ]
    net        = UCSC_NETSYNTENIC.out.net      // channel: [ val(meta), path(noClass.net.gz) ]
    axt        = UCSC_AXTSORT.out.axt          // channel: [ val(meta), path(net.axt.gz) ]
}

//
// The pair's meta map, without the tile keys of a chaining tile
//
def pairMeta(meta) {
    return meta.findAll { key, _value -> !(key in [ 'chain_tile', 'chain_tiles' ]) }
}
