#!/usr/bin/env python

###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__prog_name__ = 'compare_tk_tax.py'
__prog_desc__ = 'Compare Taxonomies from GTDB-Tk split/no_split and Official Taxonomy.'

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2021'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'p.chaumeil@uq.edu.au'
__status__ = 'Development'

import argparse
import collections
import sys
from collections import defaultdict


class Comparator(object):
    def __init__(self):
        pass

    def parse_htd(self, htd, summary_taxonomy_file, classified_high_tree_gids, gtdb_final_taxonomy,
                  genome_backbone_classification_infos, genome_metadata, red_value, genome_mapping):

        rank_table = {'d__': {'overclassified': {'genomes': [], 'rank': []},
                              'correct': {'genomes': [], 'rank': []},
                              'underclassified': {'genomes': [], 'rank': []},
                              'conflicting': {'genomes': [], 'rank': []}},
                      'p__': {'overclassified': {'genomes': [], 'rank': []},
                              'correct': {'genomes': [], 'rank': []},
                              'underclassified': {'genomes': [], 'rank': []},
                              'conflicting': {'genomes': [], 'rank': []}},
                      'c__': {'overclassified': {'genomes': [], 'rank': []},
                              'correct': {'genomes': [], 'rank': []},
                              'underclassified': {'genomes': [], 'rank': []},
                              'conflicting': {'genomes': [], 'rank': []}},
                      'o__': {'overclassified': {'genomes': [], 'rank': []},
                              'correct': {'genomes': [], 'rank': []},
                              'underclassified': {'genomes': [], 'rank': []},
                              'conflicting': {'genomes': [], 'rank': []}},
                      'f__': {'overclassified': {'genomes': [], 'rank': []},
                              'correct': {'genomes': [], 'rank': []},
                              'underclassified': {'genomes': [], 'rank': []},
                              'conflicting': {'genomes': [], 'rank': []}},
                      'g__': {'overclassified': {'genomes': [], 'rank': []},
                              'correct': {'genomes': [], 'rank': []},
                              'underclassified': {'genomes': [], 'rank': []},
                              'conflicting': {'genomes': [], 'rank': []}},
                      's__': {'overclassified': {'genomes': [], 'rank': []},
                              'correct': {'genomes': [], 'rank': []},
                              'underclassified': {'genomes': [], 'rank': []},
                              'conflicting': {'genomes': [], 'rank': []}}}

        different_cases = {}
        for item in htd:
            status = 'correct'
            rankdiff = 0
            higherrank = ''

            gid, nosplit_tax, split_rank = item
            only_high_level = 'False'
            if "_".join(gid.split("_", 2)[:2]) in classified_high_tree_gids:
                print("_".join(gid.split("_", 2)[:2]))
                only_high_level = 'True'

            for i, e in reversed(list(enumerate(split_rank.split(';')))):
                no_split_rank = nosplit_tax.split(';')
                if status != 'conflicting' and len(e) != 3 and len(no_split_rank[i]) != 3 and e == no_split_rank[i]:
                    if status != 'underclassified':
                        higherrank = e[0:3]
                    break
                elif len(e) == 3 and len(no_split_rank[i]) == 3:
                    continue

                elif status == 'conflicting' or (len(e) != 3 and len(no_split_rank[i]) != 3 and e != no_split_rank[i]):
                    status = 'conflicting'
                    if e == no_split_rank[i]:
                        break
                    higherrank = e[0:3]
                    rankdiff += 1
                elif len(e) != 3 and len(no_split_rank[i]) == 3:
                    status = 'overclassified'
                    rankdiff += 1
                    higherrank = e[0:3]
                elif len(e) == 3 and len(no_split_rank[i]) != 3:
                    if status != 'underclassified':
                        higherrank = e[0:3]
                        status = 'underclassified'
                    rankdiff += 1

            summary_taxonomy_file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t'.format(gid, split_rank, nosplit_tax,
                                                                                          gtdb_final_taxonomy.get(
                                                                                              "_".join(
                                                                                                  gid.split("_", 2)[
                                                                                                  :2]),
                                                                                              'No official taxonomy'),
                                                                                          only_high_level,
                                                                                          genome_backbone_classification_infos.get(
                                                                                              "_".join(
                                                                                                  gid.split("_", 2)[
                                                                                                  :2])),
                                                                                          genome_metadata.get(
                                                                                              "_".join(
                                                                                                  gid.split("_", 2)[
                                                                                                  :2]))[0],
                                                                                          genome_metadata.get("_".join(
                                                                                              gid.split("_", 2)[:2]))[
                                                                                              1], status, rankdiff))
            summary_taxonomy_file.write('{}\t{}'.format(red_value.get(
                gid).get('backbone'), red_value.get(gid).get('nosplit')))

            if "_".join(gid.split("_", 2)[:2]) in genome_mapping:
                summary_taxonomy_file.write('\t{}\t{}'.format(genome_mapping.get(
                    "_".join(gid.split("_", 2)[:2])).get('ani_class'),
                                                              genome_mapping.get("_".join(gid.split("_", 2)[:2])).get(
                                                                  'tree_mapped')))
            else:
                summary_taxonomy_file.write('\tN/A\tN/A')

            if 'backbone' in red_value.get(gid).keys() and 'nosplit' in red_value.get(gid).keys() and red_value.get(
                    gid).get('nosplit') != 'N/A' and red_value.get(gid).get('backbone') != 'N/A':
                summary_taxonomy_file.write('\t{}\n'.format(
                    abs(float(red_value.get(gid).get('backbone')) - float(red_value.get(gid).get('nosplit')))))
            else:
                summary_taxonomy_file.write('\tN/A\n')

            different_cases.setdefault(status, {}).setdefault(higherrank, {}).setdefault(rankdiff, []).append(
                (nosplit_tax, split_rank))

            rank_table.get(higherrank).get(
                status).get('genomes').append(item)

        # =========================================================================
        # for status, cases in different_cases.items():
        #     #print(f'Status:{status}')
        #     for frod, v in cases.items():
        #         #print(f'First rank of disagreement:{frod}')
        #         for rd,w in v.items():
        #             print(f'Status: {status}/rank of disagreement: {frod}/rank difference: {rd}')
        #             no_split_taxonomy, split_taxonomy = w[0]
        #             print("nosplit_taxonomy: {}".format(no_split_taxonomy))
        #             print("split_taxonomy: {}".format(split_taxonomy))
        #             print("------")
        #     print("\n")
        # =========================================================================

        self.create_rank_table(rank_table,gtdb_final_taxonomy)

    def create_rank_table(self, rank_table,gtdb_final_taxonomy):
        """
        Calculate statistics
        """

        status = ['correct', 'conflicting',
                  'overclassified', 'underclassified']
        print("rank\t#_genomes\tunique_taxa\t" + "\t".join(status))

        list_ranks_accross = [0] * 7
        list_ranks_accross[0] = 'Total'
        for idx_k,k in enumerate(["d__", "p__", "c__", "o__", "f__", "g__", "s__"]):
            string_result = '{}\t'.format(k)
            v = rank_table.get(k)
            nbr_genomes = 0
            rank_clade = []
            for stat in status:
                infos = v.get(stat)
                nbr_genomes += len(infos.get('genomes'))
                for gen in infos.get('genomes'):
                    rank_clade.append(gtdb_final_taxonomy.get("_".join(gen[0].split("_", 2)[:2])).split(";")[idx_k])
                    #rank_clade.append(gen[1].split(";")[idx_k])
                #print(collections.Counter(rank_clade).most_common(30))
            unique_rank_clade = set(rank_clade)
            string_result = string_result + '{}\t{}'.format(nbr_genomes, len(unique_rank_clade))


            for idx, stat in enumerate(status):
                infos = v.get(stat)
                list_ranks_accross[idx + 3] = list_ranks_accross[idx + 3] + \
                                              len(infos.get('genomes'))
                list_ranks_accross[1] = list_ranks_accross[1] + \
                                        len(infos.get('genomes'))

                # if stat.startswith('under') and k == 'o__':
                #     for x in infos.get('genomes'):
                #         if x[0] not in ['GCA_002070835.1_test','GCA_009842305.1_test','GCA_011053895.1_test','GCA_011054785.1_test','GCA_011364115.1_test','GCA_012963965.1_test','GCA_013349775.1_test','GCA_903834195.1_test','GCA_903857665.1_test','GCA_903865345.1_test','GCA_903883095.1_test','GCA_903884125.1_test','GCA_903886915.1_test','GCA_903900375.1_test','GCA_903925625.1_test','GCA_903960075.1_test','GCF_007744975.1_test']:
                #             print(x)
                #     print(len(infos.get('genomes')))
                if nbr_genomes != 0:
                    string_result = string_result + \
                                    "\t{} ({}%)".format(len(infos.get('genomes')),
                                                       round(len(infos.get('genomes')) * 100.0 / nbr_genomes, 2))
                else:
                    string_result = string_result + \
                                    "\t{} ({}%)".format(len(infos.get('genomes')), '0.0')

            print(string_result)
        print('\t'.join([str(x) for x in list_ranks_accross]))
        print("Done")

    def run_comparison(self, tk_nosplit, tk_split, high_lvl_tax, official_tax, tree_mapping,metadata_file, output_file):
        red_value = defaultdict(lambda: defaultdict(dict))

        nosplit_dict = {}
        with open(tk_nosplit, 'r') as nosplitfile:
            nosplitfile.readline()
            for line in nosplitfile:
                infos = line.strip().split('\t')
                nosplit_dict[infos[0]] = infos[1]
                red_value[infos[0]]['nosplit'] = infos[18]

        genome_backbone_classification_infos = {}
        with open(high_lvl_tax, 'r') as gltfile:
            gltfile.readline()
            for line in gltfile:
                infos = line.strip().split('\t')
                genome_backbone_classification_infos["_".join(
                    infos[0].split("_", 2)[:2])] = infos[2]
                red_value[infos[0]]['backbone'] = infos[5]

        genome_mapping = {}
        mapped_genomes = []
        with open(tree_mapping, 'r') as tmfile:
            headers = tmfile.readline().strip().split('\t')
            is_ani_classification_index = headers.index('is_ani_classification')
            species_tree_mapped_index = headers.index('species_tree_mapped')
            for line in tmfile:
                infos = line.strip().split('\t')

                genome_mapping["_".join(
                    infos[0].split("_", 2)[:2])] = {'ani_class': infos[is_ani_classification_index], 'tree_mapped': infos[species_tree_mapped_index]}
                mapped_genomes.append(infos[0])

        genome_metadata = {}
        with open(metadata_file, 'r') as ggmfile:
            headers = ggmfile.readline().strip().split('\t')
            checkm_completeness_index = headers.index('checkm_completeness')
            checkm_contamination_index = headers.index('checkm_contamination')
            for line in ggmfile:
                infos = line.strip().split('\t')
                if not infos[0].startswith('U_'):
                    genome_metadata[infos[0].replace('GB_', '').replace('RS_', '')] = (
                        infos[checkm_completeness_index], infos[checkm_contamination_index])

        split_dict = {}
        with open(tk_split, 'r') as splitfile:
            splitfile.readline()
            for line in splitfile:
                infos = line.strip().split('\t')
                split_dict[infos[0]] = infos[1]

        classified_high_tree_gids = [ "_".join(x.split("_", 2)[:2]) for x in split_dict.keys() if x not in mapped_genomes]
        #print(classified_high_tree_gids)
        #print(len(classified_high_tree_gids))

        # classified_high_tree_gids = []
        # with open(classified_only_backbone, 'r') as hightreefile:
        #     for line in hightreefile:
        #         gid = "_".join(line.split("_", 2)[:2])
        #         classified_high_tree_gids.append(gid)

        gtdb_final_taxonomy = {}
        with open(official_tax, 'r') as gtr:
            for line in gtr:
                gid, tax = line.strip().split('\t')
                gtdb_final_taxonomy[gid.replace('GB_', '').replace('RS_', '')] = tax

        #print(len(nosplit_dict), len(split_dict))

        #main_list = list(set(nosplit_dict) - set(split_dict))

        # taxonomy overclassified in split
        tois = 0

        # taxonomy underclassified in split
        tuis = 0

        # not assigned to lower tree
        natlt = 0

        # different order
        do = 0

        # higher_tree_difference
        htd = 0

        # assign_to_lower_tree_in_split
        alts = 0

        # classification different in lower tree
        ltd = 0

        countdiff = 0
        tois_list = []
        natlt_list = []
        tuis_list = []
        do_list = []
        htd_list = []
        alts_list = []
        ltd_list = []

        same_results = []

        nosplit_len = len(nosplit_dict)

        count = 0
        for k, vnosplit in nosplit_dict.items():
            count += 1
            print(f'{k} : {count}/{nosplit_len}', end='\r')
            if split_dict.get(k) != vnosplit:
                countdiff += 1
                print(split_dict.get(k))
                split_rank_list = [rank for rank in split_dict.get(k).split(';') if len(rank) > 3]
                nosplit_rank_list = [
                    rank for rank in vnosplit.split(';') if len(rank) > 3]
                # print(k, nosplit_rank_list, split_rank_list)

                if len(split_rank_list) < 4 and len(nosplit_rank_list) < 4:
                    htd += 1
                    # print("_".join(k.split("_", 2)[:2]))
                    htd_list.append(
                        (k, vnosplit, split_dict.get(k)))

                elif len(split_rank_list) > 3 and len(nosplit_rank_list) < 4:
                    alts += 1
                    alts_list.append(
                        (k, vnosplit, split_dict.get(k)))

                elif len(split_rank_list) < 4 and len(nosplit_rank_list) > 3:
                    natlt += 1
                    natlt_list.append(
                        (k, vnosplit, split_dict.get(k)))
                    # ==================================================================
                    # if nosplit_rank_list[-1][0:3] != 'o__':
                    #     print(k)
                    #     print(';'.join(nosplit_rank_list))
                    #     print(';'.join(split_rank_list))
                    # ==================================================================

                elif all(item in split_rank_list for item in nosplit_rank_list):
                    tois += 1
                    tois_list.append(
                        (k, vnosplit, split_dict.get(k)))
                elif all(item in nosplit_rank_list for item in split_rank_list):
                    tuis += 1
                    tuis_list.append(
                        (k, vnosplit, split_dict.get(k)))

                elif split_rank_list[3] != nosplit_rank_list[3]:
                    do += 1
                    do_list.append(
                        (k, vnosplit, split_dict.get(k)))
                elif all(item in nosplit_rank_list for item in
                         split_rank_list[0:4]) and nosplit_rank_list != split_rank_list:
                    ltd += 1
                    ltd_list.append(
                        (k, vnosplit, split_dict.get(k)))

                else:
                    print(k, nosplit_rank_list, split_rank_list, all(
                        item in nosplit_rank_list for item in split_rank_list[0:4]))
            else:
                same_results.append(
                    (k, vnosplit, split_dict.get(k)))

        summary_taxonomy_file = open('compare_split_nosplit.tsv', 'w')
        summary_taxonomy_file.write(
            'Genome_id\tgtdbtk_split_taxonomy\tgtdbtk_no_split_taxonomy\tgtdb_official_r202\tonly_backbone_classification\ton_terminal_branch_backbone\tcheckm_completeness\tcheckm_contamination\tstatus\trank_diff\t')
        summary_taxonomy_file.write(
            'red_value_backbone\tred_value_nosplit\tani_classification\ttree_mapping\tred_value_difference\n')

        # =========================================================================
        # higher_tree_parser = parse_htd(
        #     htd_list, summary_taxonomy_file, classified_high_tree_gids, gtdb_final_taxonomy, genome_backbone_classification_infos, genome_metadata)
        # higher_tree_parser = parse_htd(
        #     ltd_list, summary_taxonomy_file, classified_high_tree_gids, gtdb_final_taxonomy, genome_backbone_classification_infos, genome_metadata)
        # higher_tree_parser = parse_htd(
        #     natlt_list, summary_taxonomy_file, classified_high_tree_gids, gtdb_final_taxonomy, genome_backbone_classification_infos, genome_metadata)
        # higher_tree_parser = parse_htd(
        #     alts_list, summary_taxonomy_file, classified_high_tree_gids, gtdb_final_taxonomy, genome_backbone_classification_infos, genome_metadata)
        # =========================================================================
        higher_tree_parser = self.parse_htd(
            tois_list + tuis_list + do_list + alts_list + natlt_list + ltd_list + htd_list + same_results,
            summary_taxonomy_file, classified_high_tree_gids, gtdb_final_taxonomy, genome_backbone_classification_infos,
            genome_metadata, red_value, genome_mapping)
        # higher_tree_parser = parse_htd(tuis_list)
        # higher_tree_parser = parse_htd(do_list)
        summary_taxonomy_file.close()

        print('Using {} genomes'.format(len(nosplit_dict)))
        print('{} genomes have same taxonomy'.format(len(same_results)))
        print('There are {} different results:'.format(countdiff))
        print('-Both classification are higher than order level, but classification is different:')
        print('\t{}'.format(htd))
        print('-Both classification are in the same order, but classification is different:')
        print('\t{}'.format(ltd))
        print('- Genomes are not associated with a lower tree using split but are with default')
        print('\t{}'.format(natlt))
        print('- Genomes are associated with a lower tree using split but are not with default')
        print('\t{}'.format(alts))
        print('- Genomes are overclassified using split')
        print('\t{}'.format(tois))
        print('- Genomes are underclassified using split')
        print('\t{}'.format(tuis))
        print('- Genomes are classified with a different order')
        print('\t{}'.format(do))

        # =========================================================================
        # print('- {} where genome is not associated with a lower tree'.format(natlt))
        # print('(last_no_split,last_split),{}'.format(collections.Counter(natlt_list)))
        # print('- {} where genome is overclassified in Split Tree'.format(tois))
        # print('(last_no_split,last_split),{}'.format(collections.Counter(tois_list)))
        # print('- {} where genome is underclassified in Split Tree'.format(tuis))
        # print('(last_no_split,last_split),{}'.format(collections.Counter(tuis_list)))
        # print('- {} genome associated with different order.'.format(do))
        # =========================================================================


if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--tk_nosplit', required=True,
                        help='Taxonomy from GTDB-Tk with no split trees.')
    parser.add_argument('--tk_split', required=True,
                        help='Taxonomy from GTDB-Tk with split trees.')
    parser.add_argument('--high_lvl_tax', required=True,
                        help='High level taxonomy')
    parser.add_argument('--official_tax', required=True,
                        help='Final official Taxonomy.')
    parser.add_argument('--tree_mapping', required=True,
                        help='tsv file listing the Tree mapping')
    parser.add_argument('--meta_file', required=True,
                        help='Metadata file.')
    # parser.add_argument('--only_backbone', required=True,
    #                     help='lst file . Genomes on')
    parser.add_argument('--output_file', required=True,
                        help='Output file')
    args = parser.parse_args()

    try:
        gp = Comparator()
        gp.run_comparison(args.tk_nosplit, args.tk_split, args.high_lvl_tax, args.official_tax, args.tree_mapping,
                          args.meta_file,
                          args.output_file)
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise
