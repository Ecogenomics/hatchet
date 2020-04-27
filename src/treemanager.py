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
import operator

from collections import defaultdict

from biolib.common import check_file_exists, make_sure_path_exists, check_dir_exists
from biolib.seq_io import read_fasta


class TreeManager():
    def __init__(self, tree, taxonomy, rank, msa_file, domain):
        """Initialization"""
        self.logger = logging.getLogger('timestamp')

        self.tree = tree
        self.taxonomy = os.path.abspath(taxonomy)
        self.taxonomy_dict = {}
        self.phylum_dict = {}
        self.rank = rank
        self.msa_file = msa_file
        self.msa = read_fasta(msa_file)
        self.rank_order = ['d', 'p', 'c', 'o', 'f', 'g', 's']
        # we will need to select one genome from each subrank so we pick the
        # index of the subrank
        self.rank_order_index = self.rank_order.index(rank) + 1
        self.domain = 'd__Bacteria'
        if domain == 'arc':
            self.domain = 'd__Archaea'
        self.gid_per_rank = self.get_gid_per_rank(
            taxonomy, self.rank_order_index)

    def get_gid_per_rank(self, tf, roi):
        results = {}
        with open(tf) as filein:
            for line in filein:
                infos = line.strip().split('\t')
                order = infos[1].split(';')[roi]
                str_domain = infos[1].split(';')[0]
                if str_domain == self.domain:
                    results.setdefault(order, []).append(infos[0])
        return results

    def calculate_maximum_tree_size(self):
        rank_list = []
        with open(self.taxonomy) as taxf:
            for line in taxf:
                genome, taxonomy = line.strip().split('\t')
                self.taxonomy_dict[genome] = taxonomy

                ranks = taxonomy.split(';')
                if ranks[0] == self.domain:
                    rank_list.append(ranks[self.rank_order_index])
        return len(set(rank_list))

    def combine_clades(self, largest_clade_leaves, ordered_largest_clade, largest_clade_index, rank_dict_gid):
        list_gids = []
        list_ranks = []
        while largest_clade_index < len(ordered_largest_clade) and largest_clade_leaves > len(list_gids) + ordered_largest_clade[largest_clade_index][1]:
            list_gids.extend(rank_dict_gid.get(
                ordered_largest_clade[largest_clade_index][0]))
            list_ranks.append(ordered_largest_clade[largest_clade_index][0])
            largest_clade_index += 1
        return list_gids, list_ranks, largest_clade_index

    def split_tree(self, outdir):
        shell_command_file = open(os.path.join(
            outdir, 'tree_creation_commands.sh'), 'w')
        # shell_command_file.write('#!/usr/bin/env bash\n')

        list_order_file = open(os.path.join(outdir, 'order_trees.tsv'), 'w')

        make_sure_path_exists(outdir)
        mapping_file = open(os.path.join(outdir, 'tree_mapping.tsv'), 'w')
        size_maximum_tree = self.calculate_maximum_tree_size()
        self.logger.info("Maximum size is {}".format(size_maximum_tree))
        tree = dendropy.Tree.get_from_path(self.tree,
                                           schema='newick',
                                           rooting='force-rooted',
                                           preserve_underscores=True)

        # self.prepare_tree(self.tree, self.msa_file,
        #                  self.gid_per_rank, outdir, shell_command_file)

        processed_nodes = []
        tree_index = 1
        while len(tree.leaf_nodes()) > size_maximum_tree * 1.1:
            dict_nodes = {}
            print('Leaves in remaining Tree:{}.'.format(len(tree.leaf_nodes())))
            for nd in tree.levelorder_node_iter():
                if nd.is_internal() and nd not in processed_nodes:
                    if nd.label is not None and self.rank + '__' in nd.label:
                        dict_nodes[nd] = len(nd.leaf_nodes())

            largest_node = max(dict_nodes.items(),
                               key=operator.itemgetter(1))[0]

            nd_to_remove = largest_node
            processed_nodes = nd_to_remove.leaf_nodes()
            while len(nd_to_remove.parent_node.leaf_nodes()) < (size_maximum_tree * 1.1):
                nd_to_remove = nd_to_remove.parent_node
                processed_nodes = nd_to_remove.leaf_nodes()
            self.logger.info("Node to remove:{} / size of selected subtree:{} /size of parent tree:{}".format(nd_to_remove,
                                                                                                              len(nd_to_remove.leaf_nodes(
                                                                                                              )),
                                                                                                              len(nd_to_remove.parent_node.leaf_nodes())))

            def no_node_rank_filter_fn(nd): return nd.is_internal(
            ) or nd not in processed_nodes

            # Add clade_index
            self.process_tree(tree, processed_nodes,
                              tree_index, outdir, shell_command_file, list_order_file)
            tree_index += 1

            tree = tree.extract_tree(
                node_filter_fn=no_node_rank_filter_fn)

        processed_nodes = tree.leaf_nodes()

        # Add clade_index
        print("size of last tree:{}\n".format(len(tree.leaf_nodes())))

        self.process_tree(tree, processed_nodes,
                          tree_index, outdir, shell_command_file, list_order_file)

        shell_command_file.close()

    def process_tree(self, tree, processed_nodes, clade_index, outdir, shell_command_file, list_order_file):

        list_genomes = [
            nd.taxon.label for nd in processed_nodes if nd.is_leaf()]
        print("Number of genomes: {},{}".format(
            len(list_genomes), 2 * len(list_genomes)))

        set_rank_in_tree = set(
            [self.taxonomy_dict.get(gid).split(';')[self.rank_order.index(self.rank)] for gid in list_genomes])

        for rk_in_tree in set_rank_in_tree:
            list_order_file.write('{}\t{}\n'.format(rk_in_tree, clade_index))

        outgroup_id, species_out = self.select_outgroup(tree, list_genomes)
        list_genomes.append(outgroup_id)
        self.logger.info('Tree {} outgroup on {} ({})'.format(
            clade_index, species_out, outgroup_id))

        msa_file = open(os.path.join(
            outdir, '{}_msa.fa'.format(clade_index)), 'w')

        for k, v in self.msa.items():
            if k in list_genomes or k == outgroup_id:
                msa_file.write('>{}\n{}\n'.format(k, v))
        msa_file.close()

        def node_rank_filter_fn(nd): return nd.is_internal(
        ) or nd.taxon.label in list_genomes

        tree1 = tree.extract_tree(
            node_filter_fn=node_rank_filter_fn)

        package_name = 'gtdbtk.package.{}.refpkg'.format(clade_index)
        aln_file = '{}_msa.fa'.format(clade_index)
        tree_file = '{}_reference.tree'.format(clade_index)
        tree1.write_to_path(os.path.join(outdir, tree_file),
                            schema='newick',
                            suppress_rooting=True,
                            unquoted_underscores=True)

        fitted_tree_file = '{}_reference.fitted.tree'.format(clade_index)
        rooted_tree_file = '{}_reference.rooted.tree'.format(clade_index)
        stripped_tree_file = '{}_reference.stripped.tree'.format(clade_index)
        decorated_tree_file = '{}_reference.decorated.tree'.format(clade_index)

        log_fitting_file = 'log_fitting.{}.log'.format(clade_index)

        # shell_command_file.write('#!/usr/bin/env bash;')

        cmd1 = "genometreetk strip {} {}".format(
            tree_file, stripped_tree_file)
        self.purge_reload('miniconda3; module load genometreetk',
                          shell_command_file, cmd1)

        cmd2 = "FastTreeMP -wag -nome -mllen -intree {} -log {} < {} > {}".format(
            stripped_tree_file, log_fitting_file, aln_file, fitted_tree_file)
        self.purge_reload('fasttree', shell_command_file, cmd2)

        cmd3 = "genometreetk outgroup {} {} '{}' {}".format(
            fitted_tree_file, self.taxonomy, species_out, rooted_tree_file)
        self.purge_reload('genometreetk', shell_command_file, cmd3)

        cmd4 = "phylorank decorate {} {} {} --skip_rd_refine".format(
            rooted_tree_file, self.taxonomy, decorated_tree_file)
        self.purge_reload('phylorank', shell_command_file, cmd4)

        shell_command_file.write(
            "sed -i -r 's/\s+//g' {};".format(decorated_tree_file))

        cmd5 = 'taxit create -l {0} -P {0} --aln-fasta {1} --tree-stats {2} --tree-file {3}'.format(
            package_name, aln_file, log_fitting_file, decorated_tree_file)
        self.purge_reload('taxtastic/0.5.3', shell_command_file, cmd5)

        shell_command_file.write('\n')

        # sys.exit()

    def purge_reload(self, command, shell_command_file, cmd):
        shell_command_file.write('module purge;')
        shell_command_file.write('module load {};'.format(command))
        shell_command_file.write('{};'.format(cmd))

    def get_rank_for_genome(self, genome_id, rank):
        return self.taxonomy_dict.get(genome_id).split(';')[self.rank_order.index(rank)]

    def select_outgroup(self, tree, list_genomes):
        potential_outid = [nd.taxon.label for nd in tree.leaf_nodes()]
        random.seed(2)
        selected_outid = random.choice(potential_outid)
        return (selected_outid, self.get_rank_for_genome(selected_outid, 's'))
