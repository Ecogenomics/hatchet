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

__prog_name__ = 'comp_tax_new_classification.py'
__prog_desc__ = 'Compare the GTDB Taxonomony vs the GTDB-Tk Taxonomy'

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2019'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.2'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'p.chaumeil@uq.edu.au'
__status__ = 'Development'


import sys
import os
import argparse
import copy
import json


class TaxComp(object):
    def __init__(self, gtdbtk_classify_dir, tree_taxonomy, gtdb_metadata, reference_genome_file):
        self.orderrank = ["d__", "p__", "c__", "o__", "f__", "g__", "s__"]

        #self.uba_taxonomy = self.parse_gtdb_taxonomy(uba_gtdb_taxonomy)
        self.tree_taxonomy = self.parse_gtdb_taxonomy(tree_taxonomy)

        self.pplacer_confidence = {}
        self.pplacer_confidence = self.merge_two_dicts(self.parse_confidence(
            gtdbtk_classify_dir, 'ar122'), self.parse_confidence(gtdbtk_classify_dir, 'bac120'))

        # self.ar122_median_value = self.parse_red_dict(
        #     gtdbtk_classify_dir, 'ar122')
        self.bac120_median_value = self.parse_red_dict(
            gtdbtk_classify_dir, 'bac120')

        # self.ar122_gtdbtk_taxonomy = self.parse_gtdbtk_taxonomy(
        #     gtdbtk_classify_dir, 'ar122')
        self.bac120_gtdbtk_taxonomy = self.parse_gtdbtk_taxonomy(
            gtdbtk_classify_dir, 'bac120')

        # self.ar122_pplacer_taxonomy = self.parse_pplacer_taxonomy(
        #     gtdbtk_classify_dir, 'ar122')
        self.bac120_pplacer_taxonomy = self.parse_pplacer_taxonomy(
            gtdbtk_classify_dir, 'bac120')

        # self.bac120_extra_infos = self.parse_extra_infos(
        #     gtdbtk_classify_dir, 'bac120')
        # self.ar122_extra_infos = self.parse_extra_infos(
        #     gtdbtk_classify_dir, 'ar122')

        self.gtdb_metadata = self.parse_gtdb_metadata(gtdb_metadata)

        self.list_ranks, self.nospe_list_ranks = self.generate_list_ranks(
            reference_genome_file, self.tree_taxonomy)

        self.bac120_table = {'d__': {'overclassified': {'genomes': [], 'rank': []},
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
        # self.ar122_table = copy.deepcopy(self.bac120_table)

    def parse_confidence(self, gtdbtk_classify_dir, marker_id):
        result = {}
        if not os.path.isfile(os.path.join(
                gtdbtk_classify_dir, 'pplacer', 'pplacer.' + marker_id + '.json')):
            return result
        placement_hash = json.load(open(os.path.join(
            gtdbtk_classify_dir, 'pplacer', 'pplacer.' + marker_id + '.json')))
        lwr_idx = placement_hash['fields'].index('like_weight_ratio')
        for item in placement_hash['placements']:
            for node in item.get('nm'):
                id_genome = node[0]
                confs = []
                for ps in item.get('p'):
                    confs.append(str(round(ps[lwr_idx], 2)))

                id_genome = id_genome.replace('_genomic', '')
                id_genome = id_genome.replace('U_GCA_', 'GB_GCA_')

                result[id_genome] = '({})'.format(','.join(confs))
        return result

    def merge_two_dicts(self, x, y):
        z = x.copy()   # start with x's keys and values
        z.update(y)    # modifies z with y's keys and values & returns None
        return z

    def parse_extra_infos(self, gtdbtk_classify_dir, marker_id):
        result = {}
        with open(os.path.join(gtdbtk_classify_dir, 'gtdbtk.' + marker_id + '.debug_file.tsv')) as df:
            for line in df:
                info = line.rstrip().split('\t')
                id_genome = info[0].replace('_genomic', '')
                id_genome = id_genome.replace('U_GCA_', 'GB_GCA_')
                result[id_genome] = [info[2], info[3],
                                     info[4], info[5], info[6], info[7]]
        return result

    def canonical_gid(self,gid):
        """Get canonical form of NCBI genome accession.
        Example:
            G005435135 -> G005435135
            GCF_005435135.1 -> G005435135
            GCF_005435135.1_ASM543513v1_genomic -> G005435135
            RS_GCF_005435135.1 -> G005435135
            GB_GCA_005435135.1 -> G005435135
        """

        if gid.startswith('U'):
            return gid

        gid = gid.replace('RS_', '').replace('GB_', '')
        gid = gid.replace('GCA_', 'G').replace('GCF_', 'G')
        if '.' in gid:
            gid = gid[0:gid.find('.')]

        return gid

    def generate_list_ranks(self, reference_genome_file, gtdb_taxonomy):
        list_selected_genomes = []
        list_ranks = []
        with open(reference_genome_file, 'r') as rgf:
            for line in rgf:
                id = self.canonical_gid(line.rstrip())

                list_selected_genomes.append(id)

        for k, v in gtdb_taxonomy.items():
            k = self.canonical_gid(k)
            if k in list_selected_genomes:
                list_ranks.extend(v.split(";"))

        nospe_list_ranks = filter(lambda a: a != 's__', list_ranks)
        return list_ranks, nospe_list_ranks

    def parse_gtdb_metadata(self, gtdb_metadata):
        """
        Parse the metadata file to get extra information for each genome.
        """

        meta_dict = {}
        with open(gtdb_metadata, 'r') as metainfo:
            headers = metainfo.readline().split('\t')
            comp_idx = headers.index("checkm_completeness")
            conta_idx = headers.index("checkm_contamination")
            n50_contigs_index = headers.index('n50_contigs')
            contig_count_index = headers.index('contig_count')
            for line in metainfo:
                info = line.split('\t')
                k = info[0]
                if k.startswith('GB_') or k.startswith('RS_'):
                    k = k[3:]

                meta_dict[k] = {"comp": info[comp_idx], 'conta': info[conta_idx],
                                "n50": info[n50_contigs_index], "cc": info[contig_count_index]}
        return meta_dict

    def parse_pplacer_taxonomy(self, gtdbtk_classify_dir, marker_id):
        """
        Store the pplacer taxonomy for each user genome
        """
        result = {}
        with open(os.path.join(gtdbtk_classify_dir, 'classify', 'intermediate_results', 'gtdbtk.' + marker_id + '.high.classification_pplacer.tsv')) as clapf:
            for line in clapf:
                info = line.rstrip().split('\t')
                id_genome = info[0].replace('_genomic', '')
                id_genome = id_genome.replace('U_GCA_', 'GB_GCA_')
                result[id_genome] = info[1]
        return result

    def parse_gtdbtk_taxonomy(self, gtdbtk_classify_dir, marker_id):
        """
        Parse the gtdbtk summary file to store taxonomy and extra information.
        """
        results = {}
        test = []
        with open(os.path.join(gtdbtk_classify_dir, 'gtdbtk.' + marker_id + '.summary.tsv')) as sumf:
            for line in sumf:
                if line.startswith('user_genome'):
                    continue
                infos = line.rstrip().split('\t')
                id_genome = infos[0].replace('_genomic', '')
                id_genome = id_genome.replace('U_GCA_', 'GB_GCA_')
                id_genome = '_'.join(id_genome.split('_')[:2])
                test.append(infos[17])
                results[id_genome] = {
                    'tax': infos[1].replace('XXXA', ''), 'f_ref': infos[2], 'f_ani': infos[5], 'note': infos[13], 'red': infos[17]}
        return results

    def parse_red_dict(self, gtdbtk_classify_dir, marker_id):
        """
        Parse the RED dictionary.
        """

        results = {}
        results['d__'] = 0.00
        with open(os.path.join(gtdbtk_classify_dir, 'classify', 'intermediate_results', 'gtdbtk.' + marker_id + '.red_dictionary.tsv')) as dictf:
            for line in dictf:
                infos = line.strip().split('\t')
                if infos[0] == "Phylum":
                    results['p__'] = float(infos[1])
                elif infos[0] == "Class":
                    results['c__'] = float(infos[1])
                elif infos[0] == "Order":
                    results['o__'] = float(infos[1])
                elif infos[0] == "Family":
                    results['f__'] = float(infos[1])
                elif infos[0] == "Genus":
                    results['g__'] = float(infos[1])
        return results

    def parse_gtdb_taxonomy(self, gtdb_tax):
        """
        Parse the existing GTDB taxonomy 
        """
        result = {}
        duplicate = 0
        with open(gtdb_tax, 'r') as gtdbtsx:
            for line in gtdbtsx:
                info = line.rstrip().split('\t')
                k = info[0]
                if k.startswith('GB_') or k.startswith('RS_'):
                    k = k[3:]
                if k in result:
                    duplicate += 1
                result[k] = info[1].replace('XXXA', '')
        # print result.keys()
        return result

    def create_string(self, genome, infos, case, diffrank, higherdiff, gtdbtk_tax, list_ranks, gtdbtax_list, median_value, outfile):
        """
        Create a string summarising all information for each user genome.
        """
        outstring = [genome]
        outstring.append(case)
        outstring.append(str(diffrank))
        outstring.append(higherdiff)
        outstring.append(self.gtdb_metadata.get(genome).get("comp"))
        outstring.append(self.gtdb_metadata.get(genome).get("conta"))
        outstring.append(self.gtdb_metadata.get(genome).get("n50"))
        outstring.append(self.gtdb_metadata.get(genome).get("cc"))
        outstring.append(gtdbtk_tax)
        gtdbtax = self.tree_taxonomy.get(genome)
        outstring.append(";".join(self.complete_taxonomy(gtdbtax_list)))
        spelist = list(list_ranks)
        # spelist.append('s__')
        full_gtdbtax_list = [x if x in spelist else x +
                             '__x' for x in gtdbtax.split(";")]
        outstring.append(";".join(full_gtdbtax_list))
        #outstring.append(pplacer_taxonomy.get(genome))
        outstring.append('N/A')
        #======================================================================
        # EXTRA METADATA
        #======================================================================
        outstring.append(str(self.pplacer_confidence.get(genome)))

        if case != 'correct':

            if infos.get('f_ref') != 'N/A'  and float(infos.get('f_ani')) < 95.0:
                # print "ANI_issue"
                outstring.append("ANI_issue")
            # elif pplacer_taxonomy.get(genome) not in ";".join(gtdbtax_list) and infos.get('red') != 'N/A':
            #     outstring.append("pplacer_incongruency")
            elif ';'.join(gtdbtk_tax.split(";")[:len(gtdbtax_list)]) not in ";".join(gtdbtax_list) and infos.get('red') != 'N/A':
                outstring.append("pplacer_incongruency")
            elif infos.get('red') != 'N/A' and infos.get('note').endswith('topology'):
                outstring.append("pplacer_incongruency")
            elif infos.get('red') != 'N/A':
                outstring.append("RED_issue")
            else:
                outstring.append("")
        else:
            outstring.append("")

        if infos.get('red') == 'N/A':
            outstring.extend(['None', 'None', 'None', 'None', 'None', 'None',
                              'None', 'FastANI', 'None', 'None', 'None', 'None', 'None'])
        else:
            outstring.append(infos.get('red'))

            outstring.extend(
                ['None', 'None', 'None', 'None', 'None', 'None', ])
            outstring.append("RED")
            outstring.append(str(median_value.get(
                'p__') - float(infos.get('red'))))
            outstring.append(str(median_value.get(
                'c__') - float(infos.get('red'))))
            outstring.append(str(median_value.get(
                'o__') - float(infos.get('red'))))
            outstring.append(str(median_value.get(
                'f__') - float(infos.get('red'))))
            outstring.append(str(median_value.get(
                'g__') - float(infos.get('red'))))
        #======================================================================
        # END OF EXTRA METADATA
        #======================================================================

        outfile.write("{0}\n".format("\t".join(outstring)))

    def complete_taxonomy(self, tax_list):
        """
        Fill the taxonomy to have 7 ranks.
        """
        results = []
        final_results = list(self.orderrank)
        for item in tax_list:
            rank = item[0:3]
            idx = self.orderrank.index(rank)
            final_results[idx] = item
        return final_results

        tax_list.extend(self.orderrank[len(tax_list):])
        return tax_list

    def run(self, out_dir):
        """
        Main function
        """
        self.generate_taxonomy_report(
            os.path.join(
                out_dir, 'bac120_summary.tsv'), self.bac120_gtdbtk_taxonomy, self.bac120_pplacer_taxonomy, self.bac120_median_value,
            self.bac120_table, os.path.join(out_dir, 'bac120_table.tsv'))
        # self.generate_taxonomy_report(
        #     os.path.join(
        #         out_dir, 'ar122_summary.tsv'), self.ar122_gtdbtk_taxonomy, self.ar122_pplacer_taxonomy, self.ar122_median_value, self.ar122_extra_infos,
        #     self.ar122_table, os.path.join(out_dir, 'ar122_table.tsv'), self.bac120_table, os.path.join(out_dir, 'combine_table.tsv'))

    def create_rank_table(self, rk_table_f, rank_table):
        """
        Calculate statistics
        """

        status = ['correct', 'conflicting',
                  'overclassified', 'underclassified']
        print("rank\t#_genomes\t" + "\t".join(status))
        rk_table_f.write("rank\t#_genomes\t{}\n".format("\t".join(status)))

        list_ranks_accross = [0] * 4
        for idx_k, k in enumerate(self.orderrank):
            string_result = '{}\t'.format(k)
            v = rank_table.get(k)
            nbr_genomes = 0
            rank_clade = []
            for stat in status:
                infos = v.get(stat)
                nbr_genomes += len(infos.get('genomes'))
                for gen in infos.get('genomes'):
                    rank_clade.append(self.tree_taxonomy.get(gen).split(";")[idx_k])
            unique_rank_clade = set(rank_clade)
            string_result = string_result + '{}\t{}'.format(nbr_genomes, len(unique_rank_clade))

            for idx, stat in enumerate(status):
                infos = v.get(stat)
                list_ranks_accross[idx] = list_ranks_accross[idx] + \
                    len(infos.get('genomes'))
                if nbr_genomes != 0:
                    string_result = string_result + \
                        "\t{} ({}%)".format(len(infos.get('genomes')),
                                            round(len(infos.get('genomes')) * 100.0 / nbr_genomes, 2))
                else:
                    string_result = string_result + \
                        "\t{} ({}%)".format(len(infos.get('genomes')), '0.0')

            print(string_result)
            rk_table_f.write("{}\n".format(string_result))
        rk_table_f.write("Total\t{}\t{}\n".format(
            sum(list_ranks_accross), self.calculate_percentage(list_ranks_accross)))
        print(list_ranks_accross)
        print("Done")

    def calculate_percentage(self, list_ranks_accross):
        list_results = ['{} ({}%)'.format(
            x, round(x * 100.0 / sum(list_ranks_accross), 2)) for x in list_ranks_accross]
        return '\t'.join(list_results)

    def generate_taxonomy_report(self, output_file, gtdbtk_dict, pplacer_taxonomy, median_value, rank_table, rank_table_file, combine_table=None, combine_table_file=None):
        rk_table_f = open(rank_table_file, 'w')
        if combine_table:
            combine_table_f = open(combine_table_file, 'w')
        outf_summary = open(output_file, 'w')
        outf_summary.write("genome\ttk_vs_gtdb\trank_diff\thigher_rank_diff\tcompleteness\tcontamination\tn50_contigs\tcountig_count\tgtdbtk_taxonomy\tgtdb_taxonomy\tfull_gtdb_taxonomy\tpplacer_taxonomy\tpplacer_confidence\tcause_issue\tred value\thigher_rank\thigher_value\tlower_rank\tlower_value\tcase\tclosest_rank\ttool\tphylum div\tclass div\torder div\tfamily div\tgenus div\n")
        for k, v in gtdbtk_dict.items():
            gtdbtk_tax = v.get('tax')

            gtdb_tax_list = [x for x in self.tree_taxonomy.get(
                k).split(";") if x in self.list_ranks]
            gtdb_tax_list = self.complete_taxonomy(gtdb_tax_list)

            status = 'correct'
            rankdiff = 0
            higherrank = ''

            #==================================================================
            # MAIN COMPARISON
            #==================================================================
            for i, e in reversed(list(enumerate(gtdbtk_tax.split(';')))):
                if status != 'conflicting' and len(e) != 3 and len(gtdb_tax_list[i]) != 3 and e == gtdb_tax_list[i]:
                    if status != 'underclassified':
                        higherrank = e[0:3]
                    break
                elif len(e) == 3 and len(gtdb_tax_list[i]) == 3:
                    continue

                elif status == 'conflicting' or (len(e) != 3 and len(gtdb_tax_list[i]) != 3 and e != gtdb_tax_list[i]):
                    if e == gtdb_tax_list[i]:
                        break
                    if status != 'conflicting':
                        higherrank = e[0:3]

                    status = 'conflicting'
                    rankdiff += 1
                elif len(e) != 3 and len(gtdb_tax_list[i]) == 3:
                    status = 'overclassified'
                    rankdiff += 1
                    higherrank = e[0:3]
                elif len(e) == 3 and len(gtdb_tax_list[i]) != 3:
                    if status != 'underclassified':
                        higherrank = e[0:3]
                        status = 'underclassified'
                    rankdiff += 1

            self.create_string(k, v, status, rankdiff, higherrank, gtdbtk_tax,
                               self.list_ranks, gtdb_tax_list, median_value, outf_summary)

            rank_table.get(higherrank).get(
                status).get('genomes').append(k)
            if combine_table:
                combine_table.get(higherrank).get(
                    status).get('genomes').append(k)

        self.create_rank_table(rk_table_f, rank_table)
        if combine_table:
            self.create_rank_table(combine_table_f, combine_table)
            combine_table_f

        outf_summary.close()
        rk_table_f.close()


if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--gtdbtk_classify_dir',
                        help='GTDB-Tk classify directory.')
    parser.add_argument('--tree_taxonomy', help='')
    parser.add_argument('--gtdb_metadata', help='')
    parser.add_argument('--reference_genome_file', help='')

    parser.add_argument('--output_dir', help='Output directory.')

    args = parser.parse_args()

    try:
        taxcomp = TaxComp(args.gtdbtk_classify_dir, args.tree_taxonomy, args.gtdb_metadata, args.reference_genome_file)
        taxcomp.run(args.output_dir)
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise
