import logging
import os

from hatchet.biolib_lite.seq_io import read_fasta, read_seq
from hatchet.biolib_lite.common import make_sure_path_exists

from collections import defaultdict, namedtuple

logger = logging.getLogger('timestamp')


def _parse_gtdb_rep_metadata(gtdb_metadata_file, domain):
    """Parse metadata for GTDB representatives."""

    Metadata = namedtuple(
        'METADATA', 'comp cont is_mag_sag num_contigs num_ambig')

    env_assembly_types = set([
        'derived from metagenome',
        'derived from environmental sample',
        'derived from single cell'])

    gtdb_rep_metadata = {}
    with open(gtdb_metadata_file) as f:
        header = f.readline().strip().split('\t')

        gtdb_rep_idx = header.index('gtdb_representative')
        gtdb_domain_idx = header.index('gtdb_taxonomy')
        comp_idx = header.index('checkm_completeness')
        cont_idx = header.index('checkm_contamination')
        assembly_type_idx = header.index('ncbi_genome_category')
        num_contig_idx = header.index('contig_count')
        num_ambig_idx = header.index('ambiguous_bases')

        for line in f:
            tokens = line.strip().split('\t')

            is_gtdb_rep = tokens[gtdb_rep_idx].lower().startswith('t')

            if is_gtdb_rep and tokens[gtdb_domain_idx].startswith(domain):
                rid = tokens[0]

                comp = float(tokens[comp_idx])
                cont = float(tokens[cont_idx])
                assembly_type = tokens[assembly_type_idx]
                num_contigs = int(tokens[num_contig_idx])
                num_ambig = int(tokens[num_ambig_idx])

                assert assembly_type == 'none' or assembly_type in env_assembly_types

                gtdb_rep_metadata[rid] = Metadata(
                    comp,
                    cont,
                    assembly_type.lower() in env_assembly_types,
                    num_contigs,
                    num_ambig)

    return gtdb_rep_metadata


def _parse_gtdb_rep_taxonomy(gtdb_metadata_file, domain):
    """Parse GTDB taxonomy for representatives."""

    gtdb_taxonomy = {}
    with open(gtdb_metadata_file) as f:
        header = f.readline().strip().split('\t')

        gtdb_rep_idx = header.index('gtdb_representative')
        gtdb_tax_idx = header.index('gtdb_taxonomy')
        gtdb_domain_idx = header.index('gtdb_taxonomy')

        for line in f:
            tokens = line.strip().split('\t')

            is_gtdb_rep = tokens[gtdb_rep_idx].lower().startswith('t')

            if is_gtdb_rep and tokens[gtdb_domain_idx].startswith(domain):
                rid = tokens[0]
                gtdb_taxonomy[rid] = [t.strip()
                                      for t in tokens[gtdb_tax_idx].split(';')]

    return gtdb_taxonomy


def _determine_gap_percentage(gtdb_msa_file):
    """Calculate percentage of gaps in MSA."""

    gap_perc = {}
    for gid, seq in read_seq(gtdb_msa_file):
        num_gaps = sum([1 for ch in seq if ch == '-'])
        gap_perc[gid] = 100.0 * num_gaps / len(seq)

    return gap_perc

def _calculate_quality_score(gtdb_rep_metadata, gap_perc):
    """Calculate quality score for each genome.

    The quality score (QC) is based on the criteria used to determine
    the quality of genomes when updating GTDB representatives:
    https://gtdb.ecogenomic.org/faq#how-are-gtdb-species-representatives-updated-with-each-release

    Quality score = checkm_completeness
                    - 5*checkm_contamination
                    - 100*(MAG or SAG)
                    - 5*(num_contigs/100)
                    - 5*(num_undetermined_bases/10000)
                    - 5*(percentage_of_gaps_in_msa)
    """

    quality_score = {}
    for gid, cur_gap_perc in gap_perc.items():
        m = gtdb_rep_metadata[gid]

        qs = (m.comp - 5 * m.cont
              - 5 * (m.num_contigs / 100.0)
              - 5 * (m.num_ambig / 10000.0)
              - 5 * cur_gap_perc)

        if m.is_mag_sag:
            qs -= 100

        quality_score[gid] = qs

    return quality_score


def _select_high_quality_genome(rids, gtdb_rep_qual):
    """Determine highest quality genome."""

    sel_rid = None
    highest_qs = -1e6
    for rid in rids:
        if gtdb_rep_qual[rid] > highest_qs:
            highest_qs = gtdb_rep_qual[rid]
            sel_rid = rid

    return sel_rid, highest_qs

def selectbestgenomesets(metadata_file, domain, msa, rank_of_interest, output_dir):
    """Select genomes for GTDB-Tk backbone tree."""

    order_ranks = ['d', 'p', 'c', 'o', 'f', 'g', 's']

    # read metadata for each GTDB representative
    logger = logging.getLogger('timestamp')
    logger.info('Reading metadata for GTDB representatives:')
    metadata = _parse_gtdb_rep_metadata(metadata_file, domain)
    logger.info(
        f' - read metadata for {len(metadata):,} genomes')

    # read GTDB metadata for representative genomes
    logger.info('Reading GTDB taxonomy for GTDB representatives:')
    gtdb_taxonomy = _parse_gtdb_rep_taxonomy(metadata_file, domain)
    logger.info(
        f' - read taxonomy for {len(metadata):,} genomes')

    # determine percentage of gaps in MSA
    logger.info('Determining percentage of gaps in MSA:')
    gap_perc = _determine_gap_percentage(msa)
    logger.info(
        f' - determined gap percentage for {len(gap_perc):,} genomes')

    # calculate quality score for each GTDB representatives
    logger.info('Calculating quality score for GTDB representatives:')
    gtdb_rep_qual = _calculate_quality_score(metadata, gap_perc)
    logger.info(
        f' - calculated quality score for {len(gtdb_rep_qual):,} genomes')

    # get set of genomes in each taxa at the specified taxonomic rank
    logger.info('Determining genomes in each target taxon:')
    rank_labels = ['domain', 'phylum', 'class',
                   'order', 'family', 'genus', 'species']

    gids_in_taxon = defaultdict(list)
    for rid, gtdb_taxa in gtdb_taxonomy.items():
        target_taxon = gtdb_taxa[order_ranks.index(rank_of_interest)]
        gids_in_taxon[target_taxon].append(rid)

    logger.info(
        f' - identified {len(gids_in_taxon):,} {rank_labels[order_ranks.index(rank_of_interest)]}')

    # select highest quality genome from each taxon
    reps_selection_dir = os.path.join(output_dir, 'reps_selection')
    make_sure_path_exists(reps_selection_dir)
    fout = open(os.path.join(
        reps_selection_dir, f'selected_genomes_for_{rank_labels[order_ranks.index(rank_of_interest)]}.tsv'), 'w')
    fout.write('Genome ID\tTaxon\tQuality score')
    fout.write('\tCompleteness (%)\tContamination (%)')
    fout.write('\tIs MAG/SAG\tNo. contigs\tNo. ambiguous bases\tGaps (%)')
    fout.write('\tGTDB taxonomy\tCluster size\n')

    selected_genomes = {}
    for taxon, rids in gids_in_taxon.items():
        sel_rid, qs = _select_high_quality_genome(rids, gtdb_rep_qual)

        assert qs == gtdb_rep_qual[sel_rid]

        fout.write('{}\t{}\t{:.2f}\t{:.2f}\t{:.2f}\t{}\t{}\t{}\t{:.2f}\t{}\t{}\n'.format(
            sel_rid,
            taxon,
            qs,
            metadata[sel_rid].comp,
            metadata[sel_rid].cont,
            metadata[sel_rid].is_mag_sag,
            metadata[sel_rid].num_contigs,
            metadata[sel_rid].num_ambig,
            gap_perc[sel_rid],
            ';'.join(gtdb_taxonomy[sel_rid]),
            len(rids)))

        selected_genomes[taxon]={'sel_rid':sel_rid,'cluster_size':len(rids)}

    fout.close()

    return selected_genomes