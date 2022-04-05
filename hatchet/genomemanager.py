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
import logging
import random
import sys

import dendropy
import subprocess


from hatchet.biolib_lite.seq_io import read_fasta
from hatchet.tools import merge_logs, prune, remove_character,unroot

random.seed(10)


class GenomeManager():
    def __init__(self):
        """Initialization"""
        self.order_ranks = ['d', 'p', 'c', 'o', 'f', 'g', 's']
        self.ranks_dict = {'d':'domain', 'p':'phylum', 'c':'class', 'o':'order',
                            'f':'family', 'g':'genus', 's':'species'}
        self.logger = logging.getLogger('timestamp')

    def parse_taxonomy_file(self, tf, domain, roi):
        """ Parses the taxonomy file

        @param tf: str
            path to the taxonomy file
        @param domain: str
            d__Bacteria or d__Archaea
        @param roi: str
            Rank of interest in ['d', 'p', 'c', 'o', 'f', 'g', 's']

        @return: dict[str,[str]]
            Dictionary listing all genomes (values) associated with a rank of interest (key)
        """
        results = {}
        with open(tf) as filein:
            for line in filein:
                infos = line.strip().split('\t')
                taxonomy_str=infos[1].split(';')
                picked_rank = taxonomy_str[self.order_ranks.index(roi)]
                str_domain = taxonomy_str[0]
                if str_domain == domain:
                    results.setdefault(picked_rank, []).append(infos[0])

        return results

    def pick_one_genome(self, tree, msa, taxonomy_file, domain,original_log, rank_of_interest, output_dir):
        """

        @param tree: str
            Path GTDB reference tree
        @param msa: str
            Path GTDB reference MSA
        @param taxonomy_file: str
            Path GTDB Taxonomy file
        @param domain: str
            Path GTDB domain of interest
        @param original_log: str
            Path to original FastTree log file
        @param rank_of_interest: str
            rank of interest, Hatchet will pick one genome per taxa for each of this rank
        @param output_dir: str
            Path to output directory

        """
        self.logger.info(f"Picking one genome per {self.ranks_dict.get(rank_of_interest)}")
        selected_genomes = []

        if domain == 'bac':
            dom_dict = self.parse_taxonomy_file(
                taxonomy_file, 'd__Bacteria', rank_of_interest)
        else:
            dom_dict = self.parse_taxonomy_file(
                taxonomy_file, 'd__Archaea', rank_of_interest)


        # We pick randomly one genome per family
        for k, v in dom_dict.items():
            selected_genomes.append(random.choice(v))

        msa_dict = read_fasta(msa)

        selected_genome_file = open(os.path.join(output_dir, 'selected_genomes.lst'), 'w')
        selected_msa_file = open(os.path.join(output_dir, 'msa_file.fa'), 'w')
        for sg in selected_genomes:
            selected_genome_file.write('{}\n'.format(sg))
            selected_msa_file.write('>{}\n{}\n'.format(sg, msa_dict.get(sg)))
        selected_genome_file.close()
        selected_msa_file.close()

        # Taxonomy is stripped from reference tree
        subprocess.run(["genometreetk", "strip",tree,os.path.join(output_dir, 'gtdb_stripped.tree')])

        # Reference tree is pruned to only keep one genome per family
        prune(self.logger,os.path.join(output_dir, 'gtdb_stripped.tree'),
              os.path.join(output_dir, 'selected_genomes.lst'),
              os.path.join(output_dir, 'gtdb_pruned.tree'))

        # a log file is created to matches the new MSA and pruned tree
        # This log file is requiered for pplacer
        # with open(os.path.join(output_dir, 'msa_file.fa'), 'rb', 0) as in_stream, open(os.path.join(output_dir, 'fitted_tree.tree'), 'wb', 0) as out_stream:
        #     proc = subprocess.Popen(
        #         ["FastTreeMP", "-nome", "-mllen", "-intree", os.path.join(output_dir, 'gtdb_pruned.tree'), '-log',
        #          os.path.join(output_dir, 'log_fitting.log')], stdin=in_stream, stdout=out_stream)
        #     print("the commandline is {}".format(proc.args))
        #     proc.communicate()


        # we redecorate the tree
        # because the topology stays the same there should not be any polyphyletic groups
        decorated_tree = os.path.join(output_dir, 'decorated_tree.tree')
        subprocess.run(["phylorank", "decorate", os.path.join(output_dir, 'gtdb_pruned.tree'),
                        taxonomy_file,decorated_tree,"--skip_rd_refine"])

        remove_character(decorated_tree,' ')

        #we unroot the tree
        unrooted_tree = os.path.join(output_dir, 'decorated_unrooted_tree.tree')
        unrooted_undecorated_tree = os.path.join(output_dir, 'nondecorated_unrooted_tree.tree')
        if os.path.exists(unrooted_tree):
            os.remove(unrooted_tree)
        if os.path.exists(unrooted_undecorated_tree):
            os.remove(unrooted_undecorated_tree)
        unroot(decorated_tree, unrooted_tree)
        unroot(os.path.join(output_dir, 'gtdb_pruned.tree'), unrooted_undecorated_tree)



        # We are using the original log file from FastTree and use it in the pplacer package.
        # the is no fitting of the tree in the pplacer package
        # We want to keep all information from the log but we do not want to rescale the branches
        # So the idea is to replace the latest iteration from the fitting step with the original tree
        merge_logs(original_log,
                   unrooted_undecorated_tree,
                   os.path.join(output_dir, "original_merged_backbone.log"))

        # create pplacer package
        subprocess.run(["taxit","create","-l",os.path.join(output_dir,'gtdbtk_package_backbone.refpkg'),
                       "-P",os.path.join(output_dir,'gtdbtk_package_backbone.refpkg'),
                        "--aln-fasta",os.path.join(output_dir, 'msa_file.fa'),
                        "--tree-stats",os.path.join(original_log),
                        "--tree-file",unrooted_tree])


    def regenerate_red_values(self, raw_tree, pruned_tree_path, red_file, output):
        """

        @param raw_tree: str
            Path to GTDB reference tree
        @param pruned_tree_path: str
            Path to the new pruned tree
        @param red_file: str
            Path to GTDB reference RED file
        @param output:
            Path to the output file
        @return:
        """

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
        pruned_tree = dendropy.Tree.get_from_path(pruned_tree_path,
                                                  schema='newick',
                                                  rooting='force-rooted',
                                                  preserve_underscores=True)

        self._write_rd(pruned_tree, output, unpruned_tree)

    def regenerate_low_tree_red(self, split_tree_dir, raw_tree, red_file):
        """ Regenerate the red file for each species level sub tree"""
        for reffile in os.listdir(split_tree_dir):
            if reffile.endswith("_reference.tree"):
                output_file = os.path.join(split_tree_dir, 'red_value_'+reffile.replace('_reference.tree','.tsv'))
                reffile= os.path.join(split_tree_dir, reffile)
                self.regenerate_red_values(raw_tree, reffile, red_file, output_file)

    def _write_rd(self, pruned_tree, output_rd_file, unpruned_tree):
        """Write out relative divergences for each node."""

        fout = open(output_rd_file, 'w')
        for n in pruned_tree.preorder_node_iter():
            if n.is_leaf():
                fout.write('%s\t%f\n' % (n.taxon.label, 1.00))
            else:
                # get left and right taxa that define this node
                taxa = list(n.preorder_iter(lambda n: n.is_leaf()))
                # get rel_dist of this node in the original tree
                reldist_node = unpruned_tree.mrca(
                    taxon_labels=[taxa[0].taxon.label, taxa[-1].taxon.label])
                fout.write('%s|%s\t%f\n' %
                           (taxa[0].taxon.label, taxa[-1].taxon.label, reldist_node.rel_dist))

        fout.close()
