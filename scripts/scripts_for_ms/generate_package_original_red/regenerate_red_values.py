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

__prog_name__ = 'regenerate_red_values.py'
__prog_desc__ = 'Regenerate the RED values using the original branch lengths.'

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
import random
import shutil
import subprocess
import sys
from os.path import basename

import dendropy

from tools import make_sure_path_exists, prune, read_fasta, regenerate_red_values, merge_logs, symlink, yesno, unroot


class Tester(object):
    def __init__(self):
        self.rank_order = ['d', 'p', 'c', 'o', 'f', 'g', 's']
        self.sample_size = 100

        # add clade already processed to avoid re-processing
        self.processed_ranks = []


    def parse_taxonomy(self,tax_file,roi_idx):
        results = {}
        list_genomes = []
        with open(tax_file) as tf:
            for line in tf:
                infos = line.strip().split('\t')
                tax_list = infos[1].split(';')
                if tax_list[0] == 'd__Bacteria':
                    results.setdefault(tax_list[roi_idx],[]).append(infos[0])
                    list_genomes.append(infos[0])
        return results,list_genomes

    def generate_subtree(self,ref_tree,red_file,original_tk_folder,all_genomes,
                         taxonomy_number_per_rank,tax_file,msa_file,list_genomes,count,out_dir):

        msa_dict = read_fasta(msa_file)

        original_tk_folder = os.path.abspath(original_tk_folder)



        for idx,name_rank in enumerate(list_genomes):


            outdir = os.path.join(out_dir,f'{name_rank}')
            if os.path.isdir(outdir):
               print(f"We skip {name_rank}")
               continue
            make_sure_path_exists(outdir)

            subprocess.run(["genometreetk", "strip", ref_tree, os.path.join(outdir, 'gtdb_stripped.tree')])
            selected_genomes = [x for x in all_genomes if x not in taxonomy_number_per_rank.get(name_rank)]
            selected_g_file = open(os.path.join(outdir, 'selected_genomes.lst'),'w')

            submsa_file = open(os.path.join(outdir,'msa_file.fna'),'w')

            subtax_file = open(os.path.join(outdir,'sub_taxonomy.tsv'),'w')

            #parse_taxonomy
            dict_tax = {}
            with open(tax_file) as tf:
                for line in tf:
                    infos = line.strip().split('\t')
                    dict_tax[infos[0]]= infos[1]


            #genome folder to copy genomes to analyse
            genome_folder_to_analyse = os.path.join(outdir,'genomes')
            make_sure_path_exists(genome_folder_to_analyse)
            batchfile = open(os.path.join(genome_folder_to_analyse,'batchfile.tsv'),'w')

            untrimmed_msa_dict = read_fasta(os.path.join(original_tk_folder,'msa','gtdb_r95_bac120.faa'))
            untrimmed_msa_file = open(os.path.join(outdir, 'untrimmed_msa.faa'), 'w')


            genome_folder_pathfile = {}
            with open(os.path.join(original_tk_folder,'fastani','genome_paths.tsv')) as gpfile:
                for line in gpfile:
                    infos = line.strip().split(' ')
                    genome_folder_pathfile[infos[0].replace('_genomic.fna.gz','')] = infos[1]



            for gen in taxonomy_number_per_rank.get(name_rank):
                gen = gen.replace('GB_','').replace('RS_','')
                genometocopy = os.path.join(original_tk_folder,'fastani',genome_folder_pathfile.get(gen),gen+'_genomic.fna.gz')

                batchfile.write(f'{genometocopy}\t{gen}_genomic\n')
            batchfile.close()

                #shutil.copy(genometocopy,genome_folder_to_analyse)

            for selge in selected_genomes:
                selected_g_file.write(f'{selge}\n')
                submsa_file.write(f'>{selge}\n{msa_dict.get(selge)}\n')
                subtax_file.write(f'{selge}\t{dict_tax.get(selge)}\n')
                untrimmed_msa_file.write(f'>{selge}\n{untrimmed_msa_dict.get(selge)}\n')

            selected_g_file.close()
            submsa_file.close()
            subtax_file.close()
            untrimmed_msa_file.close()

            prune(os.path.join(outdir, 'gtdb_stripped.tree'),
                  os.path.join(outdir, 'selected_genomes.lst'),
                  os.path.join(outdir, 'gtdb_pruned.tree'))

            self.regenerate_red_values(ref_tree, os.path.join(outdir, 'gtdb_pruned.tree'),
                                       red_file, os.path.join(outdir, 'red_file.tsv'))

            decorated_tree = os.path.join(outdir, 'decorated_tree.tree')
            subprocess.run(["phylorank", "decorate", os.path.join(outdir, 'gtdb_pruned.tree'),
                            tax_file, decorated_tree, "--skip_rd_refine"])

            tk_package = os.path.join(outdir,'tk_package')
            no_split_pplacer_package=os.path.join(outdir,'no_split_pplacer_package')
            make_sure_path_exists(tk_package)
            make_sure_path_exists(no_split_pplacer_package)
            symlink(os.path.join(original_tk_folder,'fastani'),os.path.join(tk_package,'fastani'),True)
            symlink(os.path.join(original_tk_folder, 'markers'), os.path.join(tk_package,'markers'),True)
            symlink(os.path.join(original_tk_folder, 'masks'), os.path.join(tk_package,'masks'),True)
            symlink(os.path.join(original_tk_folder, 'metadata'), os.path.join(tk_package,'metadata'),True)
            #self.symlink(os.path.join(original_tk_folder,'mrca_red'),os.path.join(tk_package,'mrca_red'),True)
            #self.symlink(os.path.join(original_tk_folder, 'msa'), os.path.join(tk_package,'msa'),True)
            #self.symlink(os.path.join(original_tk_folder, 'pplacer'), os.path.join(tk_package,'pplacer'),True)
            symlink(os.path.join(original_tk_folder, 'radii'), os.path.join(tk_package,'radii'),True)
            #self.symlink(os.path.join(original_tk_folder, 'taxonomy'), os.path.join(tk_package,'taxonomy'),True)

            taxonomy_package = os.path.join(tk_package, 'taxonomy')
            make_sure_path_exists(taxonomy_package)
            shutil.copy(os.path.join(outdir,'sub_taxonomy.tsv'),os.path.join(taxonomy_package,'gtdb_taxonomy.tsv'))

            mrca_red_package=os.path.join(tk_package, 'mrca_red')
            make_sure_path_exists(mrca_red_package)
            symlink(os.path.join(original_tk_folder,'mrca_red','gtdbtk_r95_ar122.tsv'),os.path.join(tk_package,'mrca_red','gtdbtk_r95_ar122.tsv'),True)
            regenerate_red_values(ref_tree,decorated_tree,red_file,os.path.join(mrca_red_package,'gtdbtk_r95_bac120.tsv'))

            msa_package=os.path.join(tk_package, 'msa')
            make_sure_path_exists(msa_package)
            symlink(os.path.join(original_tk_folder,'msa','gtdb_r95_ar122.faa'),os.path.join(tk_package,'msa','gtdb_r95_ar122.faa'),True)
            shutil.copy(os.path.join(outdir, 'untrimmed_msa.faa'),os.path.join(msa_package,'gtdb_r95_bac120.faa'))

            pplacer_package=os.path.join(tk_package, 'pplacer')
            make_sure_path_exists(pplacer_package)
            symlink(os.path.join(original_tk_folder,'pplacer','gtdb_r95_ar122.refpkg'),os.path.join(tk_package,'pplacer','gtdb_r95_ar122.refpkg'),True)


            with open(os.path.join(outdir,'msa_file.fna'), 'rb', 0) as in_stream, open(
                    os.path.join(no_split_pplacer_package, 'fitted_tree.tree'), 'wb', 0) as out_stream:
                proc = subprocess.Popen(
                    ["FastTreeMP", "-nome", "-mllen", "-intree", os.path.join(outdir, 'gtdb_pruned.tree'), '-log',
                     os.path.join(no_split_pplacer_package, 'log_fitting.log')], stdin=in_stream, stdout=out_stream)
                print("the commandline is {}".format(proc.args))
                proc.communicate()

            merge_logs(os.path.join(no_split_pplacer_package, 'log_fitting.log'),
                       os.path.join(outdir, "gtdb_pruned.tree"),
                       os.path.join(no_split_pplacer_package, "log_fitting_merged.log"))

            decorated_tree = os.path.join(no_split_pplacer_package, 'decorated_tree.tree')
            subprocess.run(["phylorank", "decorate", os.path.join(outdir, 'gtdb_pruned.tree'),
                            os.path.join(outdir,'sub_taxonomy.tsv'), decorated_tree, "--skip_rd_refine"])

            # first get all lines from file
            with open(decorated_tree, 'r') as f:
                lines = f.readlines()

            # remove spaces
            lines = [line.replace(' ', '') for line in lines]

            # finally, write lines in the file
            with open(decorated_tree, 'w') as f:
                f.writelines(lines)

            # create pplacer package
            subprocess.run(["taxit", "create", "-l", os.path.join(pplacer_package, 'gtdb_r95_bac120.refpkg'),
                            "-P", os.path.join(pplacer_package, 'gtdb_r95_bac120.refpkg'),
                            "--aln-fasta", os.path.join(outdir,'msa_file.fna'),
                            "--tree-stats", os.path.join(no_split_pplacer_package, "log_fitting_merged.log"),
                            "--tree-file", decorated_tree])

            # run the split approach
            # subprocess.run(["/srv/home/uqpchaum/development/hatchet/bin/hatchet", "hatchet_wf","-d","bac",
            #                 "--ref_tree",os.path.join(outdir, 'decorated_tree.tree'),
            #                 "--msa",os.path.join(outdir,'msa_file.fna'),
            #                 '--tax',os.path.join(outdir,'sub_taxonomy.tsv'),
            #                 '-o', os.path.join(outdir,'tk_package','generate_split_tk_package'),
            #                 '--red_file',os.path.join(outdir, 'red_file.tsv')])

    def generate_backbone_red(self,reference_tree, red_file,current_dir):
        current_decorated_tree = os.path.join(current_dir,'tk_package','pplacer','gtdb_r95_bac120.refpkg','decorated_tree.tree')

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
        pruned_tree = dendropy.Tree.get_from_path(pruned_trees,
                                                  schema='newick',
                                                  rooting='force-rooted',
                                                  preserve_underscores=True)

        self._write_rd(pruned_tree, output, unpruned_tree,red_file)

    def _write_rd(self, pruned_tree, output_rd_file, unpruned_tree,red_file):
        """Write out relative divergences for each node."""

        red_infos={}
        with open(red_file) as redf:
            for line in redf:
                if '|' in line:
                    infos = line.strip().split('\t')
                    nodeone,nodetwo=infos[0].split('|')
                    red_infos[(nodeone,nodetwo)] = infos[1]

        fout = open(output_rd_file, 'w')
        lennodeider = len(list(pruned_tree.preorder_node_iter()))
        for n in tqdm(pruned_tree.preorder_node_iter(),total=lennodeider):
            if n.is_leaf():
                fout.write('%s\t%f\n' % (n.taxon.label, 1.00))
                #fout.write('%s\t%f\n' % (n.taxon.label, n.rel_dist))
            else:
                # get left and right taxa that define this node
                taxa = list(n.preorder_iter(lambda n: n.is_leaf()))
                if (taxa[0].taxon.label, taxa[-1].taxon.label) in red_infos:
                    fout.write('%s|%s\t%s\n' %
                               (taxa[0].taxon.label, taxa[-1].taxon.label, red_infos.get((taxa[0].taxon.label, taxa[-1].taxon.label))))
                else:
                    # get rel_dist of this node in the original tree
                    reldist_node = unpruned_tree.mrca(
                        taxon_labels=[taxa[0].taxon.label, taxa[-1].taxon.label])
                    fout.write('%s|%s\t%f\n' %
                               (taxa[0].taxon.label, taxa[-1].taxon.label, reldist_node.rel_dist))

        fout.close()


    def modify_red(self,ref_tree,red_file,ref_dict,dir_input):

        red_dict = ''
        with open(ref_dict, 'r') as f:
            info = f.readline().strip()
            info = info.replace('phylum','p__')
            info = info.replace('class','c__')
            info = info.replace('order','o__')
            info = info.replace('family','f__')
            info = info.replace('genus','g__')
            info = info.replace('species','s__')
            Bac120_dict = 'RED_DIST_BAC_DICT={"d__": 0.00,'+info[1:]+'\n'
        print(Bac120_dict)

        for pref in ['p__','c__','o__','f__','g__']:
            list_directory = [x for x in glob.glob(os.path.join(dir_input, f"{pref}*"))]
            if len(list_directory) == 0:
                continue
            else:
                break

        print(list_directory)
        #print(list_directory[0])

        for current_dir in list_directory:
            current_metadata = (os.path.join(current_dir,'tk_package','metadata','metadata.txt'))
            with open(current_metadata,'r') as f:
                lines = f.readlines()

            metadata_dir_path=os.path.join(current_dir,'tk_package','metadata')
            if os.path.islink(metadata_dir_path):
                os.unlink(metadata_dir_path)
            make_sure_path_exists(metadata_dir_path)
            symlink('/srv/home/uqpchaum/gtdbtk/test_for_ms/remove_and_replace_rank/original_data/release95/metadata/genome_metadata.tsv',
                    os.path.join(metadata_dir_path, 'genome_metadata.tsv'), True)
            with open(current_metadata,'w') as f:
                f.write(Bac120_dict)
                f.write(lines[1])
                f.write(lines[2])

            run_no_split_packge = True
            if run_no_split_packge:
                # replace tree in pplacer backbone package with the original one
                if os.path.exists(os.path.join(current_dir, 'gtdb_pruned.tree')):
                    os.remove(os.path.join(current_dir, 'gtdb_pruned.tree'))

                print(f'Tree pruning for {current_dir}')
                prune(os.path.join(ref_tree),
                      os.path.join(current_dir, 'selected_genomes.lst'),
                      os.path.join(current_dir, 'gtdb_pruned.tree'))
                print('Done')

                print(f'Red values for {current_dir}')
                regenerate_red_values(ref_tree, os.path.join(current_dir, 'gtdb_pruned.tree'), red_file,
                                      os.path.join(current_dir,'tk_package','mrca_red', 'gtdbtk_r95_bac120.tsv'))
                print('Done')

                decorated_tree = os.path.join(current_dir, 'decorated_tree.tree')
                if os.path.exists(decorated_tree):
                    os.remove(decorated_tree)



                # no_split_pplacer_package = os.path.join(current_dir, 'no_split_pplacer_package')
                # print(f'Merge log files for {current_dir}' )
                # merge_logs(os.path.join(no_split_pplacer_package, 'log_fitting.log'),
                #            os.path.join(current_dir, "gtdb_pruned.tree"),
                #            os.path.join(no_split_pplacer_package, "log_fitting_merged.log"))
                # print('Done')

                print(f'Decorate tree for {current_dir}')
                subprocess.run(["phylorank", "decorate", os.path.join(current_dir, 'gtdb_pruned.tree'),
                                os.path.join(current_dir,'sub_taxonomy.tsv'), decorated_tree, "--skip_rd_refine"])

                # first get all lines from file
                with open(decorated_tree, 'r') as f:
                    lines = f.readlines()

                # remove spaces
                lines = [line.replace(' ', '') for line in lines]

                # finally, write lines in the file
                with open(decorated_tree, 'w') as f:
                    f.writelines(lines)
                print('Done')


                unrooted_tree = os.path.join(current_dir, 'decorated_unrooted_tree.tree')
                unrooted_undecorated_tree = os.path.join(current_dir, 'nondecorated_unrooted_tree.tree')
                if os.path.exists(unrooted_tree):
                    os.remove(unrooted_tree)
                if os.path.exists(unrooted_undecorated_tree):
                    os.remove(unrooted_undecorated_tree)
                unroot(decorated_tree, unrooted_tree)
                unroot(os.path.join(current_dir, 'gtdb_pruned.tree'), unrooted_undecorated_tree)

                no_split_pplacer_package = os.path.join(current_dir, 'no_split_pplacer_package')
                print(f'Merge log files for {current_dir}' )
                merge_logs(os.path.join(no_split_pplacer_package, 'log_fitting.log'),
                           unrooted_undecorated_tree,
                           os.path.join(no_split_pplacer_package, "log_fitting_merged.log"))
                print('Done')



                # create pplacer package
                gtdbtk_package = os.path.join(current_dir,'tk_package','pplacer', 'gtdb_r95_bac120.refpkg')
                if os.path.exists(gtdbtk_package) and os.path.isdir(gtdbtk_package):
                    shutil.rmtree(gtdbtk_package)
                subprocess.run(["taxit", "create", "-l", os.path.join(current_dir,'tk_package','pplacer', 'gtdb_r95_bac120.refpkg'),
                                "-P", os.path.join(current_dir,'tk_package','pplacer', 'gtdb_r95_bac120.refpkg'),
                                "--aln-fasta", os.path.join(current_dir,'msa_file.fna'),
                                "--tree-stats", os.path.join(no_split_pplacer_package, "log_fitting_merged.log"),
                                "--tree-file", unrooted_tree])

            run_backbone_packge = True
            if run_backbone_packge:
                # modify the backbone package
                backbone_package = os.path.join(current_dir,'tk_package','generate_split_tk_package','high_level')

                if os.path.exists(os.path.join(backbone_package, 'gtdb_pruned.tree')):
                    os.remove(os.path.join(backbone_package, 'gtdb_pruned.tree'))

                print(f'Tree pruning for Backbone Tree {current_dir}')
                prune(os.path.join(ref_tree),
                      os.path.join(backbone_package, 'selected_genomes.lst'),
                      os.path.join(backbone_package, 'gtdb_pruned.tree'))

                print(f'Red values for Backbone Tree {current_dir}')
                regenerate_red_values(ref_tree, os.path.join(backbone_package, 'gtdb_pruned.tree'), red_file,
                                      os.path.join(current_dir,'tk_package','split','high','red','high_red_value.tsv'))
                print('Done')

                print(f'Decorate backbone tree for {current_dir}')
                decorated_tree = os.path.join(backbone_package, 'decorated_tree.tree')
                subprocess.run(["phylorank", "decorate", os.path.join(backbone_package, 'gtdb_pruned.tree'),
                                os.path.join(current_dir, 'sub_taxonomy.tsv'), decorated_tree, "--skip_rd_refine"])

                # first get all lines from file
                with open(decorated_tree, 'r') as f:
                    lines = f.readlines()

                # remove spaces
                lines = [line.replace(' ', '') for line in lines]

                # finally, write lines in the file
                with open(decorated_tree, 'w') as f:
                    f.writelines(lines)
                print('Done')

                unrooted_tree = os.path.join(backbone_package, 'decorated_unrooted_tree.tree')
                unrooted_undecorated_tree = os.path.join(backbone_package, 'nondecorated_unrooted_tree.tree')
                if os.path.exists(unrooted_tree):
                    os.remove(unrooted_tree)
                if os.path.exists(unrooted_undecorated_tree):
                    os.remove(unrooted_undecorated_tree)
                unroot(decorated_tree, unrooted_tree)
                unroot(os.path.join(backbone_package, 'gtdb_pruned.tree'), unrooted_undecorated_tree)

                print(f'Merge log files for {current_dir}')
                merge_logs(os.path.join(backbone_package, 'log_fitting.log'),
                           os.path.join(unrooted_undecorated_tree),
                           os.path.join(backbone_package, "log_fitting_merged.log"))
                print('Done')

                # create pplacer package
                gtdbtk_package = os.path.join(current_dir,'tk_package','split','high','pplacer','gtdbtk_package_high_level')
                if os.path.exists(gtdbtk_package) and os.path.isdir(gtdbtk_package):
                    shutil.rmtree(gtdbtk_package)
                subprocess.run(["taxit", "create", "-l", gtdbtk_package,
                                "-P", gtdbtk_package,
                                "--aln-fasta", os.path.join(backbone_package,'msa_file.fa'),
                                "--tree-stats", os.path.join(backbone_package, "log_fitting_merged.log"),
                                "--tree-file", unrooted_tree])

            run_species_packge = True
            if run_species_packge:
                species_package = os.path.join(current_dir, 'tk_package', 'generate_split_tk_package', 'species_level')
                list_speindex = [int(basename(x).replace('_msa.fa','')) for x in glob.glob(os.path.join(species_package, f"*_msa.fa"))]
                list_speindex.sort()
                for speindex in list_speindex:
                    selected_genomes = read_fasta(os.path.join(species_package, f"{speindex}_msa.fa"))
                    with open(os.path.join(species_package, f"{speindex}_selected_genomes.lst"), 'w') as fselected:
                        for k,v in selected_genomes.items():
                            fselected.write(f"{k}\n")
                    if os.path.exists(os.path.join(species_package, f'{speindex}_reference.tree')):
                        os.remove(os.path.join(species_package, f'{speindex}_reference.tree'))
                    print(f'Tree pruning for Species Tree {species_package}, index: {speindex}')
                    prune(os.path.join(ref_tree),
                          os.path.join(species_package, f"{speindex}_selected_genomes.lst"),
                          os.path.join(species_package, f'{speindex}_reference.tree'))
                    print('Done')

                    print(f'Red values for Species Tree {species_package}, index: {speindex}')
                    regenerate_red_values(ref_tree, os.path.join(species_package, f'{speindex}_reference.tree'), red_file,
                                          os.path.join(current_dir, 'tk_package', 'split', 'low', 'red',
                                                       f'red_value_{speindex}.tsv'))
                    print('Done')

                    decorated_tree = os.path.join(species_package, f'{speindex}_reference.decorated.tree')
                    if os.path.exists(decorated_tree):
                        os.remove(decorated_tree)

                    print(f'Decorate Species tree for {species_package}, index: {speindex}')
                    subprocess.run(["phylorank", "decorate", os.path.join(species_package, f'{speindex}_reference.tree'),
                                    os.path.join(current_dir, 'sub_taxonomy.tsv'), decorated_tree, "--skip_rd_refine"])

                    # first get all lines from file
                    with open(decorated_tree, 'r') as f:
                        lines = f.readlines()

                    # remove spaces
                    lines = [line.replace(' ', '') for line in lines]

                    # finally, write lines in the file
                    with open(decorated_tree, 'w') as f:
                        f.writelines(lines)
                    print('Done')

                    unrooted_tree = os.path.join(species_package, f'{speindex}_decorated_unrooted_tree.tree')
                    unrooted_undecorated_tree = os.path.join(species_package, f'{speindex}_nondecorated_unrooted_tree.tree')
                    if os.path.exists(unrooted_tree):
                        os.remove(unrooted_tree)
                    if os.path.exists(unrooted_undecorated_tree):
                        os.remove(unrooted_undecorated_tree)
                    unroot(decorated_tree, unrooted_tree)
                    unroot(os.path.join(species_package, f'{speindex}_reference.tree'), unrooted_undecorated_tree)

                    print(f'Merge log files for Species Tree {species_package}, index: {speindex}')
                    merge_logs(os.path.join(species_package, f'log_fitting.{speindex}.log'),
                               unrooted_undecorated_tree,
                               os.path.join(species_package, f"log_fitting_merged.{speindex}.log"))
                    print('Done')


                    # create pplacer package
                    gtdbtk_package = os.path.join(current_dir, 'tk_package', 'split', 'low', 'pplacer',
                                                  f'gtdbtk.package.{speindex}.refpkg')
                    if os.path.exists(gtdbtk_package) and os.path.isdir(gtdbtk_package):
                        shutil.rmtree(gtdbtk_package)
                    subprocess.run(["taxit", "create", "-l", gtdbtk_package,
                                    "-P", gtdbtk_package,
                                    "--aln-fasta", os.path.join(species_package, f"{speindex}_msa.fa"),
                                    "--tree-stats", os.path.join(species_package, f"log_fitting_merged.{speindex}.log"),
                                    "--tree-file", unrooted_undecorated_tree])





if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--ref_tree', required=True,
                        help='reference tree.')
    parser.add_argument('--red_file', required=True,
                        help='RED file')
    parser.add_argument('--red_dict', required=True,
                        help='Dictionary of average RED values')
    parser.add_argument('--dir_input', required=True,
                        help='Directory listing a set of clades starting with the known prefixes.')

    args = parser.parse_args()

    try:
        tt = Tester()
        tt.modify_red(args.ref_tree,args.red_file,args.red_dict,args.dir_input)
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise