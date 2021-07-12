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

from hatchet.biolib_lite.common import make_sure_path_exists
from hatchet.biolib_lite.seq_io import read_fasta


class TreeManager():
    def __init__(self, tree, taxonomy, rank, msa_file, domain):
        """Initialization"""
        self.logger = logging.getLogger('timestamp')

        self.tree = tree
        self.taxonomy = os.path.abspath(taxonomy)
        self.taxonomy_dict = {}
        self.phylum_dict = self.populate_phylum_dict()
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

    def populate_phylum_dict(self):
        results = {}

        temp_tree = dendropy.Tree.get_from_path(self.tree,
                                           schema='newick',
                                           rooting='force-rooted',
                                           preserve_underscores=True)
        list_reps = [nd.taxon.label for nd in temp_tree.leaf_nodes()]
        with open(self.taxonomy) as taxf:
            for line in taxf:
                infos = line.strip().split('\t')
                if infos[0] not in list_reps:
                    continue
                domain,phylum,*_ = infos[1].split(';')
                if domain == 'd__Bacteria':
                    results.setdefault(phylum,[]).append(infos[0])
        return results



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

        shell_command_file.write("""

alias hatchet='/srv/home/uqpchaum/development/hatchet/bin/hatchet'
        
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/opt/centos7/sw/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/opt/centos7/sw/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/opt/centos7/sw/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/opt/centos7/sw/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
\n""")


        #list_order_file = open(os.path.join(outdir, 'tree_mapping.tsv'), 'w')

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
            self.logger.info("size of selected subtree:{} /size of parent tree:{}".format(len(nd_to_remove.leaf_nodes(
                                                                                                              )),
                                                                                                              len(nd_to_remove.parent_node.leaf_nodes())))

            def no_node_rank_filter_fn(nd): return nd.is_internal(
            ) or nd not in processed_nodes

            # Add clade_index
            self.process_tree(processed_nodes,
                              tree_index, outdir, shell_command_file, mapping_file)
            tree_index += 1

            tree = tree.extract_tree(
                node_filter_fn=no_node_rank_filter_fn)

        processed_nodes = tree.leaf_nodes()

        # Add clade_index
        print("size of last tree:{}\n".format(len(tree.leaf_nodes())))

        self.process_tree(processed_nodes,
                          tree_index, outdir, shell_command_file, mapping_file)

        shell_command_file.close()

    def select_remaining_phylum(self, set_phylum_in_tree):
        results = {}
        random.seed(2)
        for k,v in self.phylum_dict.items():
            if k in set_phylum_in_tree:
                continue
            results[random.choice(v)] = len(v)
        self.logger.info(f'Remaining Phylum = {len(results)}, Phylum in Tree {len(set_phylum_in_tree)} {set_phylum_in_tree}')
        return results

    def process_tree(self, processed_nodes, clade_index, outdir, shell_command_file, mapping_file):

        list_genomes = [
            nd.taxon.label for nd in processed_nodes if nd.is_leaf()]
        # print("Number of genomes: {},{}".format(
        #     len(list_genomes), 2 * len(list_genomes)))

        set_rank_in_tree = set(
            [self.taxonomy_dict.get(gid).split(';')[self.rank_order.index(self.rank)] for gid in list_genomes])
        set_phylum_in_tree = set(
            [self.taxonomy_dict.get(gid).split(';')[1] for gid in list_genomes])

        for rk_in_tree in set_rank_in_tree:
            mapping_file.write('{}\t{}\n'.format(rk_in_tree, clade_index))

        # Pick one genome for each Phylum not present in the process nodes
        external_phylum_genomes = self.select_remaining_phylum(set_phylum_in_tree)

        outgroup_id, species_out = self.select_outgroup(external_phylum_genomes)
        self.logger.info('Tree {} outgroup on {} ({})'.format(
            clade_index, species_out, outgroup_id))

        msa_file = open(os.path.join(
            outdir, '{}_msa.fa'.format(clade_index)), 'w')

        count_msa = 0
        taxonomy_msa = {}
        for k, v in self.msa.items():
            if k in list_genomes or k in external_phylum_genomes:
                taxonomy_msa[k] = self.taxonomy_dict.get(k)
                count_msa += 1
                msa_file.write('>{}\n{}\n'.format(k, v))
        msa_file.close()

        def node_rank_filter_fn(nd): return nd.is_internal(
        ) or nd.taxon.label in list_genomes+list(external_phylum_genomes.keys())

        sub_tree = dendropy.Tree.get_from_path(self.tree,
                                           schema='newick',
                                           rooting='force-rooted',
                                           preserve_underscores=True)

        tree1 = sub_tree.extract_tree(
            node_filter_fn=node_rank_filter_fn)

        nber_of_leaves = len(tree1.leaf_nodes())
        taxonomy_leaves = {}
        for nd1 in tree1.leaf_nodes():
            taxonomy_leaves[nd1.taxon.label] = self.taxonomy_dict.get(nd1.taxon.label)



        self.logger.info(f'Number of seq in MSA: {count_msa} / Number of leaves in Tree : {nber_of_leaves}')
        if count_msa != nber_of_leaves:
            self.logger.error("There is a difference between the number of leaves and numbers of sequences in MSA")
            keys1 = taxonomy_msa.keys()
            keys2 = taxonomy_leaves.keys()
            diff_keys = keys1-keys2
            for df in diff_keys:
                self.logger.error(f"{df}\t{self.taxonomy_dict.get(df)}")

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
        log_fitting_merged_file = 'log_fitting_merged.{}.log'.format(clade_index)

        # shell_command_file.write('#!/usr/bin/env bash;')

        cmd1 = "genometreetk strip {} {}".format(
            tree_file, stripped_tree_file)
        self.purge_reload(
                          shell_command_file, cmd1,'conda activate genometreetk-0.1.6;')

        cmd2 = "FastTreeMP -wag -nome -mllen -intree {} -log {} < {} > {}".format(
            stripped_tree_file, log_fitting_file, aln_file, fitted_tree_file)
        self.purge_reload(shell_command_file, cmd2)



        cmd_hatchet = f"hatchet merge_logs -i {log_fitting_file} --pruned_tree {stripped_tree_file} -o {log_fitting_merged_file}"
        self.purge_reload(shell_command_file, cmd_hatchet,'conda deactivate;conda activate gtdbtk-dev;')

        cmd3 = "genometreetk outgroup {} {} '{}' {}".format(
            stripped_tree_file, self.taxonomy, species_out, rooted_tree_file)
        self.purge_reload(shell_command_file, cmd3,'conda deactivate;conda activate genometreetk-0.1.6;')

        cmd4 = "phylorank decorate {} {} {} --skip_rd_refine".format(
            rooted_tree_file, self.taxonomy, decorated_tree_file)
        self.purge_reload(shell_command_file, cmd4,'conda activate phylorank-0.1.9;')

        shell_command_file.write(
            "sed -i -r 's/\s+//g' {};".format(decorated_tree_file))

        cmd5 = 'taxit create -l {0} -P {0} --aln-fasta {1} --tree-stats {2} --tree-file {3}'.format(
            package_name, aln_file, log_fitting_merged_file, decorated_tree_file)
        self.purge_reload(shell_command_file, cmd5,'conda activate taxtastic-0.9.0;')

        shell_command_file.write('\n')

        # sys.exit()

    def purge_reload(self, shell_command_file, cmd,conda_activate = None):
        if conda_activate is not None:
            shell_command_file.write(conda_activate)
        shell_command_file.write('{};'.format(cmd))

    def get_rank_for_genome(self, genome_id, rank):
        return self.taxonomy_dict.get(genome_id).split(';')[self.rank_order.index(rank)]

    # def select_outgroup(self, tree, list_genomes):
    #     potential_outid = [nd.taxon.label for nd in tree.leaf_nodes()]
    #     random.seed(2)
    #     selected_outid = random.choice(potential_outid)
    #     return (selected_outid, self.get_rank_for_genome(selected_outid, 's'))

    def select_outgroup(self, list_genomes):
        sorted_size_phylum = {k: v for k, v in sorted(list_genomes.items(), key=lambda item: item[1])}
        selected_outid = list(sorted_size_phylum.keys())[0]
        return (selected_outid, self.get_rank_for_genome(selected_outid, 's'))
