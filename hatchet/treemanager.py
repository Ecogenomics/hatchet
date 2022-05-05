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
import subprocess
import sys
import csv
import re
import logging
import random
import dendropy
import operator

from collections import defaultdict, Counter

from hatchet.biolib_lite.common import make_sure_path_exists
from hatchet.biolib_lite.seq_io import read_fasta
from hatchet.tools import merge_logs, remove_character, unroot


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
        self.rank_order_index = self.rank_order.index(rank)
        self.rank_family_index = self.rank_order.index(rank) + 1


        self.domain = 'd__Bacteria'
        if domain == 'arc':
            self.domain = 'd__Archaea'

    def populate_phylum_dict(self):
        """
            Generates a dictionary listing all representatives genomes (values) per phylum (key)

        @return: dict
            dictionary listing all representatives genomes (values) per phylum (key)
        """
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
        """
            Count how many orders are in the tree
            This number represent the size of the order level tree , i.e. backbone tree


        @return: int
            Return the size of the backbone tree
        """
        rank_list = []
        order_list = []
        with open(self.taxonomy) as taxf:
            for line in taxf:
                genome, taxonomy = line.strip().split('\t')
                self.taxonomy_dict[genome] = taxonomy

                ranks = taxonomy.split(';')
                if ranks[0] == self.domain:
                    rank_list.append(ranks[self.rank_family_index])
                    order_list.append(ranks[self.rank_order_index])

        most_common_order,len_order = Counter(order_list).most_common(1)[0]


        print(most_common_order, len_order)

        return max(len(set(rank_list)),len_order)

    def combine_clades(self, largest_clade_leaves, ordered_largest_clade, largest_clade_index, rank_dict_gid):
        list_gids = []
        list_ranks = []
        while largest_clade_index < len(ordered_largest_clade) and largest_clade_leaves > len(list_gids) + ordered_largest_clade[largest_clade_index][1]:
            list_gids.extend(rank_dict_gid.get(
                ordered_largest_clade[largest_clade_index][0]))
            list_ranks.append(ordered_largest_clade[largest_clade_index][0])
            largest_clade_index += 1
        return list_gids, list_ranks, largest_clade_index

    def split_tree(self,original_log, outdir):

        make_sure_path_exists(outdir)
        mapping_file = open(os.path.join(outdir, 'tree_mapping.tsv'), 'w')
        size_maximum_tree = self.calculate_maximum_tree_size()
        self.logger.info("Maximum size is {}".format(size_maximum_tree))
        tree = dendropy.Tree.get_from_path(self.tree,
                                           schema='newick',
                                           rooting='force-rooted',
                                           preserve_underscores=True)

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
            self.process_tree(processed_nodes,original_log,
                              tree_index, outdir, mapping_file)
            tree_index += 1

            tree = tree.extract_tree(
                node_filter_fn=no_node_rank_filter_fn)

        processed_nodes = tree.leaf_nodes()

        # Add clade_index
        print("size of last tree:{}\n".format(len(tree.leaf_nodes())))

        self.process_tree(processed_nodes,original_log,
                          tree_index, outdir, mapping_file)


    def select_remaining_phylum(self, set_phylum_in_tree):
        results = {}
        random.seed(2)
        for k,v in self.phylum_dict.items():
            if k in set_phylum_in_tree:
                continue
            results[random.choice(v)] = len(v)
        self.logger.info(f'Remaining Phylum = {len(results)}, Phylum in Tree {len(set_phylum_in_tree)} {set_phylum_in_tree}')
        return results

    def process_tree(self, processed_nodes,original_log, clade_index, outdir, mapping_file):

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
        self.logger.info(f'This number represents the number of genomes in the order(s) of interest '
                         f'({len(list_genomes)} genomes and one genome per phylum not classifiying '
                         f'order of interest ({len(external_phylum_genomes)} genomes).')

        if count_msa != nber_of_leaves:
            self.logger.error("There is a difference between the number of leaves and numbers of sequences in MSA")
            keys1 = taxonomy_msa.keys()
            keys2 = taxonomy_leaves.keys()
            diff_keys = keys1-keys2
            for df in diff_keys:
                self.logger.error(f"{df}\t{self.taxonomy_dict.get(df)}")
            sys.exit(-1)

        package_name = os.path.join(outdir,'gtdbtk.package.{}.refpkg'.format(clade_index))
        aln_file = os.path.join(outdir,'{}_msa.fa'.format(clade_index))
        tree_file = os.path.join(outdir,'{}_reference.tree'.format(clade_index))
        tree1.write_to_path(tree_file,
                            schema='newick',
                            suppress_rooting=True,
                            unquoted_underscores=True)

        fitted_tree_file = os.path.join(outdir,'{}_reference.fitted.tree'.format(clade_index))
        rooted_tree_file = os.path.join(outdir,'{}_reference.rooted.tree'.format(clade_index))
        stripped_tree_file = os.path.join(outdir,'{}_reference.stripped.tree'.format(clade_index))
        decorated_tree_file = os.path.join(outdir,'{}_reference.decorated.tree'.format(clade_index))

        #log_fitting_file = os.path.join(outdir,'log_fitting.{}.log'.format(clade_index))
        log_fitting_merged_file = os.path.join(outdir,'original_merged.{}.log'.format(clade_index))


        # We strip the taxonomy from the tree
        subprocess.run(["genometreetk", "strip", tree_file, stripped_tree_file])


        # a log file is created to matches the new MSA and pruned tree
        # This log file is requiered for pplacer
        # with open(aln_file, 'rb', 0) as in_stream, open(fitted_tree_file, 'wb', 0) as out_stream:
        #     proc = subprocess.Popen(
        #         ["FastTreeMP", "-nome", "-mllen", "-intree", stripped_tree_file, '-log',log_fitting_file], stdin=in_stream, stdout=out_stream)
        #     print("the commandline is {}".format(proc.args))
        #     proc.communicate()



        subprocess.run(["genometreetk","outgroup",stripped_tree_file,
                        self.taxonomy,species_out,rooted_tree_file])


        # we redecorate the tree
        # because the topology stays the same there should not be any polyphyletic groups
        subprocess.run(["phylorank","decorate",rooted_tree_file,
                        self.taxonomy,decorated_tree_file,"--skip_rd_refine"])

        remove_character(decorated_tree_file,' ')

        #we unroot the tree
        unrooted_tree = os.path.join(outdir, '{}_decorated_unrooted_tree.tree'.format(clade_index))
        unrooted_undecorated_tree = os.path.join(outdir, '{}_nondecorated_unrooted_tree.tree'.format(clade_index))
        if os.path.exists(unrooted_tree):
            os.remove(unrooted_tree)
        if os.path.exists(unrooted_undecorated_tree):
            os.remove(unrooted_undecorated_tree)
        unroot(decorated_tree_file, unrooted_tree)
        unroot(rooted_tree_file, unrooted_undecorated_tree)


        # This is a step when we modify the log fitting file from pplacer
        # We want to keep all information from the log but we do not want to rescale the branches
        # So the idea is to replace the latest iteration from the fitting step with the original tree
        # merge_logs(original_log,
        #            unrooted_undecorated_tree,
        #            log_fitting_merged_file)

        # create pplacer package
        subprocess.run(["taxit","create","-l",package_name,
                       "-P",package_name,
                        "--aln-fasta",aln_file,
                        "--tree-stats",original_log,
                        "--tree-file",unrooted_tree])


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
