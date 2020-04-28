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

import os
import sys
import csv
import re
import logging
import random
import dendropy
import glob
from collections import defaultdict


from biolib.common import check_file_exists, make_sure_path_exists, check_dir_exists
from biolib.seq_io import read_fasta

from src.tools import purge_reload

class GenomeManager():
    def __init__(self):
        """Initialization"""
        self.order_ranks = ['d', 'p', 'c', 'o', 'f', 'g', 's']
        self.logger = logging.getLogger('timestamp')

    def parse_taxonomy_file(self, tf, domain, roi):
        results = {}
        with open(tf) as filein:
            for line in filein:
                infos = line.strip().split('\t')
                picked_rank = infos[1].split(';')[self.order_ranks.index(roi)]
                str_domain = infos[1].split(';')[0]
                if str_domain == domain:
                    results.setdefault(picked_rank, []).append(infos[0])

        return results

    def pick_one_genome(self, tree, msa, taxonomy_file, domain, rank_of_interest, output_dir):
        selected_genomes = []
        dom_dict = {}

        if domain == 'bac':
            dom_dict = self.parse_taxonomy_file(
                taxonomy_file, 'd__Bacteria', rank_of_interest)
        else:
            dom_dict = self.parse_taxonomy_file(
                taxonomy_file, 'd__Archaea', rank_of_interest)

        for k, v in dom_dict.items():
            selected_genomes.append(random.choice(v))

        msa_dict = read_fasta(msa)

        selected_genome_file = open(os.path.join(
            output_dir, 'selected_genomes.lst'), 'w')
        selected_msa_file = open(os.path.join(output_dir, 'msa_file.fa'), 'w')
        for sg in selected_genomes:
            selected_genome_file.write('{}\n'.format(sg))
            selected_msa_file.write('>{}\n{}\n'.format(sg, msa_dict.get(sg)))
        selected_genome_file.close()
        selected_msa_file.close()

        # Because module calls are not working properly in Python we generate a
        # sh script to run after

        # We pruned the tree
        sh_file = open(os.path.join(output_dir, 'pick_one_genome.sh'), 'w')
        cmd = 'genometreetk strip {} {}'.format(
            tree, os.path.join(output_dir, 'gtdb_stripped.tree'))
        sh_file.write('module purge\n')
        sh_file.write('module load genometreetk\n')
        sh_file.write(cmd + '\n')

        cmd = 'genetreetk prune --tree {} -t {} --output {}'.format(os.path.join(output_dir, 'gtdb_stripped.tree'),
                                                                    os.path.join(
                                                                        output_dir, 'selected_genomes.lst'),
                                                                    os.path.join(output_dir, 'gtdb_pruned.tree'))
        sh_file.write('module purge\n')
        sh_file.write('module load genetreetk\n')
        sh_file.write(cmd + '\n')

        cmd = 'FastTreeMP -nome -mllen -intree {} -log {} < {} > {}'.format(os.path.join(output_dir, 'gtdb_pruned.tree'),
                                                                            os.path.join(
                                                                                output_dir, 'log_fitting.log'),
                                                                            os.path.join(
                                                                                output_dir, 'msa_file.fa'),
                                                                            os.path.join(output_dir, 'fitted_tree.tree'))


        sh_file.write(cmd + '\n')


        decorated_tree = os.path.join(output_dir, 'decorated_tree.tree')
        cmd4 = "phylorank decorate {} {} {} --skip_rd_refine".format(
            os.path.join(output_dir, 'fitted_tree.tree'),
            taxonomy_file, decorated_tree)
        purge_reload('phylorank', sh_file, cmd4)

        sh_file.write(
            "sed -i -r 's/\s+//g' {};\n".format(decorated_tree))


        # create pplacer package
        sh_file.write('module purge\n')
        sh_file.write('module load taxtastic/0.5.3\n')
        cmd5 = 'taxit create -l {0} -P {0} --aln-fasta {1} --tree-stats {2} --tree-file {3}\n'.format(
            os.path.join(output_dir,'gtdbtk_package_high_level'),
            os.path.join(output_dir, 'msa_file.fa'),
            os.path.join(output_dir, 'log_fitting.log'), decorated_tree)
        sh_file.write(cmd5 + '\n')

    def regenerate_red_values(self, raw_tree, pruned_trees, red_file, output):
        unpruned_tree = dendropy.Tree.get_from_path(raw_tree,
                                                    schema='newick',
                                                    rooting='force-rooted',
                                                    preserve_underscores=True)

        # create map from leave labels to tree nodes
        leaf_node_map = {}
        for leaf in unpruned_tree.leaf_node_iter():
            leaf_node_map[leaf.taxon.label] = leaf

        # parse RED file and associate reference RED value to reference node in
        # the tree
        reference_nodes = set()
        with open(red_file) as rf:
            for line in rf:
                label_ids, red_value = line.strip().split('\t')
                labels = label_ids.split('|')

                if len(labels) == 2:
                    if all(elem in leaf_node_map for elem in labels):
                        taxa = [leaf_node_map[label].taxon for label in labels]
                        node = unpruned_tree.mrca(taxa=taxa)
                        node.rel_dist = float(red_value)
                        reference_nodes.add(node)
                elif len(labels) == 1:
                    if labels[0] in leaf_node_map:
                        node = leaf_node_map[labels[0]]
                        node.rel_dist = float(red_value)
                        reference_nodes.add(node)
        print(pruned_trees)
        pruned_tree = dendropy.Tree.get_from_path(pruned_trees,
                                                  schema='newick',
                                                  rooting='force-rooted',
                                                  preserve_underscores=True)

        self._write_rd(pruned_tree, output, unpruned_tree)

    def regenerate_low_tree_red(self, split_tree_dir, raw_tree, red_file):
        for reffile in os.listdir(split_tree_dir):
            if reffile.endswith("_reference.tree"):
                output_file = os.path.join(split_tree_dir, reffile.replace(
                    '_reference.tree', '_red_value.txt'))
                reffile= os.path.join(split_tree_dir, reffile)
                self.regenerate_red_values(raw_tree, reffile, red_file, output_file)

    def _write_rd(self, pruned_tree, output_rd_file, unpruned_tree):
        """Write out relative divergences for each node."""

        fout = open(output_rd_file, 'w')
        for n in pruned_tree.preorder_node_iter():
            if n.is_leaf():
                fout.write('%s\t%f\n' % (n.taxon.label, 1.00))
                #fout.write('%s\t%f\n' % (n.taxon.label, n.rel_dist))
            else:
                # get left and right taxa that define this node
                taxa = list(n.preorder_iter(lambda n: n.is_leaf()))
                # get rel_dist of this node in the original tree
                reldist_node = unpruned_tree.mrca(
                    taxon_labels=[taxa[0].taxon.label, taxa[-1].taxon.label])
                fout.write('%s|%s\t%f\n' %
                           (taxa[0].taxon.label, taxa[-1].taxon.label, reldist_node.rel_dist))

        fout.close()
