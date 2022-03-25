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

__prog_name__ = 'pick_one_clade.py'
__prog_desc__ = 'Remove genome from the reference tree and try to place them again to see if they get the same taxonomy.'

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2021'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'p.chaumeil@uq.edu.au'
__status__ = 'Development'

import argparse
import os
import random
import shutil
import subprocess
import sys


from tools import make_sure_path_exists, prune, read_fasta, regenerate_red_values, merge_logs,symlink,yesno



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


    def run_tests(self,ref_tree,red_file,msa_file,genome_folder,tax_file,
                  roi,out_dir):

        random.seed(10)

        roi_idx = self.rank_order.index(roi)

        taxonomy_number_per_rank,all_genomes = self.parse_taxonomy(tax_file,roi_idx)
        genome_per_name= {}
        for k,v in taxonomy_number_per_rank.items():
            if len(v) == 1:
                key='1'
            else:
                key='many'
            genome_per_name.setdefault(key,[]).append(k)
            genome_per_name.setdefault('all',[]).append(k)

        # We removed processed ranks
        list_all_g = [x for x in (genome_per_name.get('all')) if x not in self.processed_ranks]
        list_all_g = (random.sample(list_all_g, self.sample_size))
        print(list_all_g)

        ans = yesno("Happy with this list?")
        if ans is False:
            sys.exit()


        self.generate_subtree(ref_tree,red_file,genome_folder,all_genomes,taxonomy_number_per_rank,tax_file,msa_file,list_all_g,'all',out_dir)


if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--ref_tree', required=True,
                        help='reference tree.')
    parser.add_argument('--msa_file', required=True,
                        help='msa file')
    parser.add_argument('--tax_file', required=True,
                        help='Taxonomy file')
    parser.add_argument('--red_file', required=True,
                        help='RED file')
    parser.add_argument('--roi','--rank_of_interest', required=True, choices=['d', 'p', 'c', 'o', 'f', 'g', 's'] ,
                        help='Rank to remove and generate package ')
    parser.add_argument('--original_tk_folder', required=True,
                        help='Genome folder where all genome sequences are located.')
    parser.add_argument('--out_dir', required=True,
                        help='Output directory')

    args = parser.parse_args()

    try:
        tt = Tester()
        tt.run_tests(args.ref_tree,args.red_file, args.msa_file, args.original_tk_folder,
                     args.tax_file,args.roi, args.out_dir)
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise