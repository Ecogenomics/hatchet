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

__prog_name__ = 'check_pplacer_placement.py'
__prog_desc__ = 'Check if genomes is place on the same branch between reference tree, and gtdb-tk trees.'

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2021'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'p.chaumeil@uq.edu.au'
__status__ = 'Development'

import argparse
import glob
import os
import sys
from shutil import copy

import dendropy

from tools import make_sure_path_exists, prune, read_fasta, regenerate_red_values, merge_logs

sys.setrecursionlimit(15000)

class Tester(object):
    def __init__(self):
        self.rank_order = ['d', 'p', 'c', 'o', 'f', 'g', 's']

    def list_ref_children_nodes(self,child_node,list_genomes_ref):
        return [subleaf for subleaf in child_node.leaf_iter()
                                        if subleaf.taxon.label.replace("'", '') not in list_genomes_ref]

    def list_ref_children_in_backbone(self,child_node, list_selected_genomes):
        return [subleaf for subleaf in child_node.leaf_iter()
                                        if subleaf.taxon.label.replace("'", '') in list_selected_genomes]


    def measure_split_placement_backbone(self,backbone_tree,reference_tree,list_genomes,list_genomes_ref,roi):
        # measure nosplit
        map_pruned_tree = []
        results = {}
        map_pruned_tree_reverse = {}
        for idx,nd in enumerate(backbone_tree.preorder_node_iter()):
            print(f'enumerate backbone_tree {idx}     ', end='\r')
            if nd.parent_node is None:
                child_nodes = nd.child_nodes()
                childnd_one=self.list_ref_children_nodes(child_nodes[0],list_genomes)
                childnd_two=self.list_ref_children_nodes(child_nodes[1],list_genomes)
                map_pruned_tree.append(nd)

            elif nd.is_internal():
                child_nodes = nd.child_nodes()
                childnd_one=self.list_ref_children_nodes(child_nodes[0],list_genomes)
                childnd_two=self.list_ref_children_nodes(child_nodes[1],list_genomes)

                if len(childnd_one) > 0 and len(childnd_two) > 0:
                    # children_list = childnd_one+childnd_two
                    # temp_parent_node = nd.parent_node
                    # while len(self.list_ref_children_nodes(temp_parent_node,list_genomes)) == len(children_list):
                    #     temp_parent_node = temp_parent_node.parent_node

                    map_pruned_tree.append(nd)

        list_selected_genomes = [tempnd.taxon.label for tempnd in self.list_ref_children_nodes(backbone_tree.seed_node,list_genomes)]

        map_reference_tree ={}
        cc=0
        leaf_parent_node = None
        child_one_two = []
        for idx,leaf in enumerate(reference_tree.leaf_node_iter()):
            if leaf.taxon.label in list_genomes_ref:
                cc +=1
                print(f'enumerate reference_tree {cc}/{len(list_genomes)}                   ', end='\r')
                if leaf_parent_node is None:
                    leaf_parent_node = leaf.parent_node
                    child_nodes = leaf_parent_node.child_nodes()
                    childnd_one = self.list_ref_children_in_backbone(child_nodes[0],list_selected_genomes)
                    childnd_two = self.list_ref_children_in_backbone(child_nodes[1],list_selected_genomes)
                    while len(childnd_one) == 0 or len(childnd_two) == 0:
                        leaf_parent_node = leaf_parent_node.parent_node
                        child_nodes = leaf_parent_node.child_nodes()
                        childnd_one=self.list_ref_children_in_backbone(child_nodes[0],list_selected_genomes)
                        childnd_two=self.list_ref_children_in_backbone(child_nodes[1],list_selected_genomes)
                        child_one_two = [tempnd.taxon.label for tempnd in self.list_ref_children_in_backbone(leaf_parent_node,list_selected_genomes)]

                map_reference_tree[leaf.taxon.label] = child_one_two

        same_placement = 0
        cc = 0
        pruned_parent_node = None
        pruned_pr_start = None
        for leaf in backbone_tree.leaf_node_iter():
            # print(leaf.taxon.label)

            if leaf.taxon.label in list_genomes:
                cc += 1
                print(f'backbone_tree leaf_node_iter {cc}/{len(list_genomes)}          ',end='\r')
                temp_label = leaf.taxon.label.replace('_genomic.fna', '').replace('GCA_', 'GB_GCA_').replace('GCF_',
                                                                                                             'RS_GCF_')
                pruned_parent_node = leaf.parent_node
                # if pruned_parent_node is None:
                while pruned_parent_node not in map_pruned_tree:
                    pruned_parent_node = pruned_parent_node.parent_node

                # we have the parent node on the pruned tree
                # we need to find the real parent node based on the reference tree
                mrca = backbone_tree.mrca(taxon_labels=map_reference_tree.get(temp_label))
                mrca_start = backbone_tree.mrca(taxon_labels=map_reference_tree.get(temp_label))

                # if pruned_pr_start is None:
                #     pruned_pr_start = pruned_parent_node.clone(2)
                if mrca == pruned_parent_node:
                    same_placement += 1
                    results[temp_label] = f'{roi}\t{temp_label}\t{len(list_genomes_ref)}\t0'

                else:
                    distance = 0
                    ref_parent_nodes = []
                    ref_parent_nodes_length = {}
                    while mrca is not None and mrca != pruned_parent_node:
                        if len(ref_parent_nodes) == 0:
                            ref_parent_nodes.append(mrca)
                            ref_parent_nodes_length[mrca]=len(self.list_ref_children_nodes(mrca, list_genomes))
                        elif mrca in map_pruned_tree:
                            ref_parent_nodes.append(mrca)
                        mrca = mrca.parent_node
                    if mrca == pruned_parent_node:
                        distance = len(ref_parent_nodes)
                    pruned_parent_nodes = []
                    while pruned_parent_node is not None and mrca_start != pruned_parent_node:
                        # we dont add parent node if one branch points to reference genomes
                        if len(pruned_parent_nodes) == 0:
                            pruned_parent_nodes.append(pruned_parent_node)
                        elif (len(pruned_parent_nodes) > 0 and
                              len(self.list_ref_children_nodes(pruned_parent_node, list_genomes)) != len(
                                    self.list_ref_children_nodes(pruned_parent_nodes[-1], list_genomes))):
                            pruned_parent_nodes.append(pruned_parent_node)
                        pruned_parent_node = pruned_parent_node.parent_node
                    if pruned_parent_node == mrca_start:
                        distance = len(pruned_parent_nodes)
                    if distance == 0:
                        parent_nodes = pruned_parent_nodes + ref_parent_nodes
                        distance = len([x for x in parent_nodes if parent_nodes.count(x) == 1]) + 1
                    results[temp_label] = f'{roi}\t{temp_label}\t{len(list_genomes_ref)}\t{distance}'


        return results

    def measure_split_placement(self,species_level_trees,reference_tree,list_genomes,list_genomes_ref,roi):
        results = {}
        list_genomes_processed = list_genomes.copy()
        for species_level_tree_file,sub_list_genomes in species_level_trees.items():
            species_level_tree = dendropy.Tree.get_from_path(species_level_tree_file,
                                                      schema='newick',
                                                      rooting='force-rooted',
                                                      preserve_underscores=True)

            sub_list_genomes_ref = [x.replace('_genomic.fna', '').replace('GCA_', 'GB_GCA_').replace('GCF_', 'RS_GCF_') for
                                x in sub_list_genomes]
            list_reference_genomes = []
            map_pruned_tree = []


            for idx, nd in enumerate(species_level_tree.preorder_node_iter()):
                if nd.parent_node is None:
                    child_nodes = nd.child_nodes()
                    childnd_one = self.list_ref_children_nodes(child_nodes[0], list_genomes)
                    childnd_two = self.list_ref_children_nodes(child_nodes[1], list_genomes)
                    map_pruned_tree.append(nd)
                    list_reference_genomes=[x.taxon.label for x in childnd_one + childnd_two]

                elif nd.is_internal():
                    child_nodes = nd.child_nodes()
                    childnd_one = self.list_ref_children_nodes(child_nodes[0], list_genomes)
                    childnd_two = self.list_ref_children_nodes(child_nodes[1], list_genomes)

                    if len(childnd_one) > 0 and len(childnd_two) > 0:
                        children_list = childnd_one + childnd_two
                        temp_parent_node = nd.parent_node
                        while len(self.list_ref_children_nodes(temp_parent_node, list_genomes)) == len(children_list):
                            temp_parent_node = temp_parent_node.parent_node

                        map_pruned_tree.append(nd)


            same_placement = 0
            cc = 0
            species_parent_node_dict = {}
            for leaf in species_level_tree.leaf_node_iter():
                if leaf.taxon.label in sub_list_genomes:
                    cc += 1
                    print(f'{cc}/{len(sub_list_genomes)}', end='\r')
                    temp_label = leaf.taxon.label.replace('_genomic.fna', '').replace('GCA_', 'GB_GCA_').replace('GCF_',
                                                                                                                 'RS_GCF_')
                    species_parent_node = leaf.parent_node
                    while species_parent_node not in map_pruned_tree:
                        species_parent_node = species_parent_node.parent_node
                    species_parent_node_dict[temp_label]=species_parent_node


            map_reference_tree = {}
            species_parent_node_start = None
            for idx, leaf in enumerate(reference_tree.leaf_node_iter()):
                if leaf.taxon.label in sub_list_genomes_ref:
                    temp_label = leaf.taxon.label
                    leaf_parent_node = leaf.parent_node
                    while len(self.list_ref_children_nodes(leaf_parent_node, sub_list_genomes_ref)) == 0:
                        leaf_parent_node = leaf_parent_node.parent_node

                    list_taxon_start = [x.taxon.label for x in self.list_ref_children_nodes(leaf_parent_node, sub_list_genomes_ref)]
                    if not all(elem in list_reference_genomes for elem in list_taxon_start):
                        print(f'{temp_label} not in the right species tree')
                        results[temp_label] = f'{roi}\t{temp_label}\t{len(list_genomes_ref)}\t1000'
                    else:

                        # the current leaf_parent_node has one branch with reference genomes but not the others so it means it will be removed in the pruned tree
                        # we go one level up
                        leaf_parent_node = leaf_parent_node.parent_node
                        list_taxon = [x.taxon.label for x in leaf_parent_node.leaf_nodes() if
                                      x.taxon.label in list_reference_genomes]
                        while len(list_taxon_start) == len(list_taxon):
                            leaf_parent_node = leaf_parent_node.parent_node
                            list_taxon= [x.taxon.label for x in leaf_parent_node.leaf_nodes() if
                             x.taxon.label in list_reference_genomes]

                        mrca = species_level_tree.mrca(taxon_labels=list_taxon)
                        mrca_start = species_level_tree.mrca(taxon_labels=list_taxon)

                        species_parent_node = species_parent_node_dict.get(leaf.taxon.label)

                        # if species_parent_node_start is None:
                        #     species_parent_node_start = species_parent_node.clone(2)

                        if mrca == species_parent_node:
                            same_placement +=1
                            results[temp_label] = f'{roi}\t{temp_label}\t{len(list_genomes_ref)}\t0'

                            #outputfile.write(f'{roi}\t{temp_label}\t{len(list_genomes_ref)}\t0\n')

                        else:
                            distance = 0
                            ref_parent_nodes = []
                            while mrca is not None and mrca != species_parent_node:
                                if len(ref_parent_nodes) == 0:
                                    ref_parent_nodes.append(mrca)
                                elif (len(ref_parent_nodes) > 0 and
                                        len(self.list_ref_children_nodes(mrca,list_genomes_ref)) != len(self.list_ref_children_nodes(ref_parent_nodes[-1],list_genomes_ref))):
                                    ref_parent_nodes.append(mrca)
                                mrca = mrca.parent_node
                            if mrca == species_parent_node:
                                distance = len(ref_parent_nodes)

                            leaf_parent_nodes =[]
                            while species_parent_node is not None and mrca_start != species_parent_node:
                                # we dont add parent node if one branch points to reference genomes
                                if len(leaf_parent_nodes) == 0 :
                                    leaf_parent_nodes.append(species_parent_node)
                                elif (len(leaf_parent_nodes) > 0 and
                                        len(self.list_ref_children_nodes(species_parent_node,list_genomes)) != len(self.list_ref_children_nodes(leaf_parent_nodes[-1],list_genomes))):
                                    leaf_parent_nodes.append(species_parent_node)
                                species_parent_node = species_parent_node.parent_node
                            if species_parent_node == mrca_start:
                                distance = len(leaf_parent_nodes)
                            if distance == 0:
                                parent_nodes = leaf_parent_nodes+ref_parent_nodes
                                distance = len([x for x in parent_nodes if parent_nodes.count(x) == 1]) +1
                            results[temp_label]= f'{roi}\t{temp_label}\t{len(list_genomes_ref)}\t{distance}'

        for x in list_genomes:
            temp_label = x.replace('_genomic.fna', '').replace('GCA_', 'GB_GCA_').replace('GCF_','RS_GCF_')
            if temp_label not in results:
                results[temp_label] = f'{roi}\t{temp_label}\t{len(list_genomes_ref)}\t1000'
        return results

    def measure_nosplit_placement(self,pruned_tree,reference_tree,list_genomes,list_genomes_ref,roi,results):
        # measure nosplit
        map_pruned_tree = []
        map_pruned_tree_reverse = {}
        temp_parent_node = None
        for idx,nd in enumerate(pruned_tree.preorder_node_iter()):
            print(f'enumerate pruned_tree {idx}',end='\r')
            if nd.parent_node is None:
                child_nodes = nd.child_nodes()
                childnd_one=self.list_ref_children_nodes(child_nodes[0],list_genomes)
                childnd_two=self.list_ref_children_nodes(child_nodes[1],list_genomes)
                map_pruned_tree.append(nd)

            elif nd.is_internal():
                child_nodes = nd.child_nodes()
                childnd_one=self.list_ref_children_nodes(child_nodes[0],list_genomes)
                childnd_two=self.list_ref_children_nodes(child_nodes[1],list_genomes)

                if len(childnd_one) > 0 and len(childnd_two) > 0:
                    # children_list = childnd_one+childnd_two
                    # if temp_parent_node is None:
                    #     temp_parent_node = nd.parent_node
                    #     while len(self.list_ref_children_nodes(temp_parent_node,list_genomes)) == len(children_list):
                    #         temp_parent_node = temp_parent_node.parent_node

                    map_pruned_tree.append(nd)

        map_reference_tree ={}
        cc=0
        leaf_parent_node = None
        temp_list = []
        for idx,leaf in enumerate(reference_tree.leaf_node_iter()):
            if leaf.taxon.label in list_genomes_ref:
                cc+=1
                print(f'enumerate reference_tree {cc}/{len(list_genomes_ref)}              ',end='\r')
                if leaf_parent_node is None:
                    leaf_parent_node = leaf.parent_node
                    temp_list=self.list_ref_children_nodes(leaf_parent_node,list_genomes_ref)
                    while len(temp_list) == 0:
                        leaf_parent_node = leaf_parent_node.parent_node
                        temp_list = self.list_ref_children_nodes(leaf_parent_node, list_genomes_ref)
                    #the current leaf_parent_node has one branch with reference genomes but not the others so it means it will be removed in the pruned tree
                    # we go one level up
                    if leaf_parent_node.parent_node is not None:
                        leaf_parent_node = leaf_parent_node.parent_node
                        temp_list = self.list_ref_children_nodes(leaf_parent_node, list_genomes_ref)

                map_reference_tree[leaf.taxon.label] = [x.taxon.label for x in temp_list]

        same_placement=0
        cc = 0
        pruned_parent_node = None
        pruned_pr_start = None
        for leaf in pruned_tree.leaf_node_iter():

            if leaf.taxon.label in list_genomes:
                cc+=1
                print(f'enumrate pruned_tree {cc}/{len(list_genomes)}                  ',end='\r')
                temp_label = leaf.taxon.label.replace('_genomic.fna', '').replace('GCA_', 'GB_GCA_').replace('GCF_', 'RS_GCF_')
                pruned_parent_node = leaf.parent_node
                #if pruned_parent_node is None:
                while pruned_parent_node not in map_pruned_tree:
                    pruned_parent_node = pruned_parent_node.parent_node

            # we have the parent node on the pruned tree
            # we need to find the real parent node based on the reference tree

                mrca = pruned_tree.mrca(taxon_labels=map_reference_tree.get(temp_label))
                mrca_start = mrca

                # pruned_tree.mrca(taxon_labels=map_reference_tree.get(temp_label))

                # if pruned_pr_start is None:
                #     pruned_pr_start = pruned_parent_node.clone(2)

                #print(leaf.taxon.label,mrca,pruned_parent_node)

                if mrca == pruned_parent_node:
                    same_placement +=1
                    results[temp_label]=results.get(temp_label)+(f'\t0\n')

                else:
                    distance = 0
                    ref_parent_nodes = []
                    while mrca is not None and mrca != pruned_parent_node:
                        if len(ref_parent_nodes) == 0:
                            ref_parent_nodes.append(mrca)
                        elif len(ref_parent_nodes) > 0 and mrca in map_pruned_tree:
                            ref_parent_nodes.append(mrca)
                        mrca = mrca.parent_node
                    if mrca == pruned_parent_node:
                        distance = len(ref_parent_nodes)
                    pruned_parent_nodes =[]
                    while pruned_parent_node is not None and mrca_start != pruned_parent_node:
                        # we dont add parent node if one branch points to reference genomes
                        if len(pruned_parent_nodes) == 0 :
                            pruned_parent_nodes.append(pruned_parent_node)
                        elif len(pruned_parent_nodes) > 0 and pruned_parent_node in map_pruned_tree:
                            pruned_parent_nodes.append(pruned_parent_node)
                        pruned_parent_node = pruned_parent_node.parent_node
                    if pruned_parent_node == mrca_start:
                        distance = len(pruned_parent_nodes)
                    if distance == 0:
                        parent_nodes = pruned_parent_nodes+ref_parent_nodes
                        distance = len([x for x in parent_nodes if parent_nodes.count(x) == 1]) +1
                    results[temp_label] = results.get(temp_label) + f'\t{distance}\n'


        return results

    def run(self, ref_tree,roi,nosplit_tree,split_tree_directory,placed_genomes_dir,output_file):
        outputfile = open(output_file,'w')
        outputfile.write('Rank\tGenome\tgenome pool\tdistance split\tdistance nosplit\n')

        print(f'for {roi}')

        #List genomes to check
        list_genomes = []

        for file in glob.glob(os.path.join(placed_genomes_dir,"*.fna.gz")):
            list_genomes.append((os.path.basename(file).replace('.gz','')))

        high_level_tree = dendropy.Tree.get_from_path(os.path.join(split_tree_directory,'gtdbtk.high.bac120.classify.tree'),
                                                    schema='newick',
                                                    rooting='force-rooted',
                                                    preserve_underscores=True)


        species_level_tree = {}

        for tree in glob.glob(os.path.join(split_tree_directory, "gtdbtk.bac120.classify.tree.*")):
            with open(tree) as t:
                tree_string = t.read()
                for lg in list_genomes:
                    if lg in tree_string:
                        species_level_tree.setdefault(tree,[]).append(lg)

        pruned_tree = dendropy.Tree.get_from_path(nosplit_tree,
                                                    schema='newick',
                                                    rooting='force-rooted',
                                                    preserve_underscores=True)

        reference_tree = dendropy.Tree.get_from_path(ref_tree,
                                                     schema='newick',
                                                     rooting='force-rooted',
                                                     preserve_underscores=True)



        list_genomes_ref = [x.replace('_genomic.fna', '').replace('GCA_', 'GB_GCA_').replace('GCF_', 'RS_GCF_') for x in
                            list_genomes]



        # Measure Split Placement
        print("Begin Split Parsing.....")
        if roi[0] in ['f','g']:
            results = self.measure_split_placement(species_level_tree,reference_tree,list_genomes,list_genomes_ref,roi)
        elif roi[0] in ['p', 'c', 'o']:
            results = self.measure_split_placement_backbone(high_level_tree,reference_tree,list_genomes,list_genomes_ref,roi)
        print("End Split Parsing.....")

        # Measure No Split Placement
        print("Begin NO Split Parsing.....")
        results = self.measure_nosplit_placement(pruned_tree,reference_tree,list_genomes,list_genomes_ref,roi,results)
        print("End NO Split Parsing.....")



        for k,v in results.items():
            outputfile.write(v)
        outputfile.close()


        outputfile.close()


if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-r','--ref_tree', required=True,
                        help='reference tree.')
    parser.add_argument('-n','--nosplit_tree', required=True,
                        help='GTDB-Tk No split tree.')
    parser.add_argument('-s','--split_tree_directory', required=True,
                        help='GTDB-Tk split tree directory.')
    parser.add_argument('-g','--placed_genomes_dir', required=True,
                        help='directory of "user" genomes replaced in GTDB-Tk')
    parser.add_argument('--rank', required=True,
                        help='Rank of Interest')
    parser.add_argument('-o','--output_file', required=True,
                        help='Output file')

    args = parser.parse_args()

    try:
        tt = Tester()
        tt.run(args.ref_tree,args.rank,args.nosplit_tree,args.split_tree_directory,args.placed_genomes_dir,args.output_file)
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise



    # def run(self, ref_tree,roi,classify_tree,placed_genomes_dir,output_file):
    #
    #     print(f'for {roi}')
    #
    #     #List genomes to check
    #     list_genomes = []
    #     for file in glob.glob(os.path.join(placed_genomes_dir,"*.fna.gz")):
    #         list_genomes.append((os.path.basename(file).replace('.gz','')))
    #
    #     pruned_tree = dendropy.Tree.get_from_path(classify_tree,
    #                                                 schema='newick',
    #                                                 rooting='force-rooted',
    #                                                 preserve_underscores=True)
    #
    #     reference_tree = dendropy.Tree.get_from_path(ref_tree,
    #                                                  schema='newick',
    #                                                  rooting='force-rooted',
    #                                                  preserve_underscores=True)
    #
    #     list_genomes_ref = [x.replace('_genomic.fna', '').replace('GCA_', 'GB_GCA_').replace('GCF_', 'RS_GCF_') for x in
    #                         list_genomes]
    #
    #     map_reference_tree = {}
    #     for nd in reference_tree.postorder_node_iter():
    #         if nd
    #         map_reference_tree[nd]=[childnd for childnd  in nd.postorder_iter() if childnd.taxon is not None and childnd.taxon.label not in list_genomes_ref]
    #
    #     for individual_genomes in list_genomes:
    #         genome_pruned_node_child = []
    #         genome_pruned_node_parent = []
    #         for leaf in pruned_tree.leaf_node_iter():
    #             #print(leaf.taxon.label)
    #             if leaf.taxon.label == individual_genomes:
    #                 parent_node = leaf.parent_node
    #                 while len([subleaf for subleaf in parent_node.leaf_iter()
    #                                     if subleaf.taxon.label.replace("'", '') not in list_genomes]) < 1:
    #                     parent_node = parent_node.parent_node
    #
    #                 genome_pruned_node_child = [subleaf for subleaf in parent_node.leaf_iter()
    #                                     if subleaf.taxon.label.replace("'", '') not in list_genomes]
    #
    #                 while len([subleaf for subleaf in parent_node.leaf_iter()
    #                                     if subleaf.taxon.label.replace("'", '') not in list_genomes]) == len(genome_pruned_node_child):
    #                     parent_node = parent_node.parent_node
    #                 genome_pruned_node_parent = [subleaf for subleaf in parent_node.leaf_iter()
    #                                     if subleaf.taxon.label.replace("'", '') not in list_genomes]
    #
    #
    #         genome_ref_node_child = []
    #         genome_ref_node_parent = []
    #         for refleaf in reference_tree.leaf_node_iter():
    #             #print(leaf.taxon.label)
    #
    #             if refleaf.taxon.label == individual_genomes.replace('_genomic.fna','').replace('GCA_','GB_GCA_').replace('GCF_','RS_GCF_'):
    #                 refparent_node = refleaf.parent_node
    #                 while len([subleaf for subleaf in refparent_node.leaf_iter()
    #                                     if subleaf.taxon.label.replace("'", '') not in list_genomes_ref]) < 1:
    #                     refparent_node = refparent_node.parent_node
    #
    #                 genome_ref_node_child = [subleaf for subleaf in refparent_node.leaf_iter()
    #                                     if subleaf.taxon.label.replace("'", '') not in list_genomes_ref]
    #
    #                 while len([subleaf for subleaf in refparent_node.leaf_iter()
    #                                     if subleaf.taxon.label.replace("'", '') not in list_genomes_ref]) == len(genome_ref_node_child):
    #                     refparent_node = refparent_node.parent_node
    #                 genome_ref_node_parent = [subleaf for subleaf in refparent_node.leaf_iter()
    #                                     if subleaf.taxon.label.replace("'", '') not in list_genomes_ref]
    #
    #
    #         if len(genome_ref_node_child) != len(genome_pruned_node_child) or len(genome_ref_node_parent) != len(genome_pruned_node_parent):
    #
    #             print(individual_genomes.replace('_genomic.fna','').replace('GCA_','GB_GCA_').replace('GCF_','RS_GCF_'))
    #             print(len(genome_pruned_node_child), len(genome_pruned_node_parent))
    #             print(len(genome_ref_node_child),len(genome_ref_node_parent))
    #
    #             mrca_start = reference_tree.mrca(taxon_labels=[genome_pruned_node_parent[0].taxon.label,
    #                                                      genome_pruned_node_parent[-1].taxon.label])
    #             mrca = reference_tree.mrca(taxon_labels=[genome_pruned_node_parent[0].taxon.label,
    #                                                      genome_pruned_node_parent[-1].taxon.label])
    #
    #             one_list_branches = []
    #             print(mrca)
    #             while refparent_node != mrca and mrca is not None:
    #                 one_list_branches.append(mrca)
    #                 mrca=mrca.parent_node
    #             print(mrca)
    #
    #             two_list_branches = []
    #             print(refparent_node,mrca_start,refparent_node != mrca_start)
    #             while refparent_node not in [mrca_start,mrca] and refparent_node not in one_list_branches and refparent_node is not None:
    #                 two_list_branches.append(refparent_node)
    #                 refparent_node=refparent_node.parent_node
    #             print(refparent_node)
    #
    #
    #             print(one_list_branches)
    #             print(two_list_branches)
    #
    #             if refparent_node in [mrca_start,mrca]:
    #                 print(len(two_list_branches) + 1)
    #             elif refparent_node :
    #                 print(len(two_list_branches)+one_list_branches.index(refparent_node)+1)
    #             print('###')
