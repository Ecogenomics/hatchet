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

__prog_name__ = 'pick_one_genome_per_order.py'
__prog_desc__ = 'Pick one genome per order.'

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2019'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'p.chaumeil@uq.edu.au'
__status__ = 'Development'

import os
import sys
import argparse
import random

from biolib.common import make_sure_path_exists
from biolib.seq_io import read_fasta


class genomePicker(object):
    def __init__(self, uba_mapping_file, taxonomy_file):
        """Initialization."""
        self.uba_mappings = self.parse_uba_mapping_file(uba_mapping_file)
        self.taxonomy_bact_dict = self.parse_taxonomy_file(
            taxonomy_file, 'd__Bacteria')
        self.taxonomy_arch_dict = self.parse_taxonomy_file(
            taxonomy_file, 'd__Archaea')

    def parse_uba_mapping_file(self, lgf):
        results = []
        with open(lgf) as filein:
            for line in filein:
                infos = line.strip().split('\t')
                results.append(infos[0])
                results.append(infos[1])
                if infos[2].startswith('GCA'):
                    results.append(infos[2])
                    results.append('GB_' + infos[2])
                results.append(infos[0])
                results.append(line.strip())
        return results

    def parse_taxonomy_file(self, tf, domain):
        results = {}
        with open(tf) as filein:
            for line in filein:
                infos = line.strip().split('\t')
                if infos[0] == 'GB_GCA_002254385.1':
                    continue
                family = infos[1].split(';')[4]
                str_domain = infos[1].split(';')[0]
                if str_domain == domain:
                    if family in results:
                        results.get(family).append(infos[0])
                    else:
                        results[family] = [infos[0]]
        return results

    def prepare_tree(self, tree, msa_file, dom_dict, output_dir):
        selected_genomes = []
        for k, v in dom_dict.iteritems():
            selected_genomes.append(random.choice(v))

        msa_dict = read_fasta(msa_file)

        selected_genome_file = open(os.path.join(
            output_dir, 'selected_genomes.lst'), 'w')
        selected_msa_file = open(os.path.join(output_dir, 'msa_file.fa'), 'w')
        for sg in selected_genomes:
            selected_genome_file.write('{}\n'.format(sg))
            selected_msa_file.write('>{}\n{}\n'.format(sg, msa_dict.get(sg)))
        selected_genome_file.close()
        selected_msa_file.close()

        # We pruned the tree
        cmd = 'genometreetk strip {} {}'.format(
            tree, os.path.join(output_dir, 'gtdb_stripped.tree'))
        print('module purge\n')
        print('module load genometreetk')
        print(cmd)
        # os.system(cmd)

        cmd = 'genetreetk prune --tree {} -t {} --output {}'.format(os.path.join(output_dir, 'gtdb_stripped.tree'),
                                                                    os.path.join(
                                                                        output_dir, 'selected_genomes.lst'),
                                                                    os.path.join(output_dir, 'gtdb_pruned.tree'))
        print('module purge\n')
        print('module load genetreetk')
        print(cmd)
        # os.system(cmd)

        cmd = 'FastTreeMP -nome -mllen -intree {} -log {} < {} > {}'.format(os.path.join(output_dir, 'gtdb_pruned.tree'),
                                                                            os.path.join(
                                                                                output_dir, 'log_fitting.log'),
                                                                            os.path.join(
                                                                                output_dir, 'msa_file.fa'),
                                                                            os.path.join(output_dir, 'fitted_tree.tree'))
        print(cmd)
        # os.system(cmd)

    def run(self, bac_tree, bac_msa, arc_tree, arc_msa, output_dir):
        make_sure_path_exists(output_dir)
        make_sure_path_exists(os.path.join(output_dir, 'bacteria'))
        make_sure_path_exists(os.path.join(output_dir, 'archaea'))
        self.prepare_tree(bac_tree, bac_msa, self.taxonomy_bact_dict,
                          os.path.join(output_dir, 'bacteria'))
        self.prepare_tree(arc_tree, arc_msa, self.taxonomy_arch_dict,
                          os.path.join(output_dir, 'archaea'))
        return True


if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--bac_tree', required=True,
                        help='GTDB-Tk Tree .')
    parser.add_argument('--bac_msa', required=True,
                        help='GTDB-Tk MSA file.')
    parser.add_argument('--arc_tree', required=True,
                        help='GTDB-Tk Tree .')
    parser.add_argument('--arc_msa', required=True,
                        help='GTDB-Tk MSA file.')
    parser.add_argument('--taxonomy', required=True,
                        help='GTDB-Tk Taxonomy file.')
    parser.add_argument('--uba_mappings', required=True,
                        help='')
    parser.add_argument('--output_dir', required=True,
                        help='')
    args = parser.parse_args()

    try:
        gp = genomePicker(args.uba_mappings, args.taxonomy)
        gp.run(args.bac_tree, args.bac_msa, args.arc_tree,
               args.arc_msa, args.output_dir)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
