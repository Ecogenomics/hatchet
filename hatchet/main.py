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
import logging

from hatchet.biolib_lite.common import make_sure_path_exists
from hatchet.common import package_results
from hatchet.tools import merge_logs, unroot

from hatchet.genomemanager import GenomeManager
from hatchet.treemanager import TreeManager


class OptionsParser():
    def __init__(self):
        """Initialization"""

        self.logger = logging.getLogger('timestamp')

    def pick_genomes(self, options):
        make_sure_path_exists(options.out_dir)
        p = GenomeManager()
        p.pick_one_genome(options.ref_tree,options.metadata,
                          options.msa, options.tax,
                          options.domain,options.original_log,
                          options.rank_of_interest, options.out_dir)

    def regenerate_red_values(self, options):
        p = GenomeManager()
        p.regenerate_red_values(
            options.ref_tree, options.pruned_tree, options.red_file, options.out_file)

    def regenerate_low_tree_red(self, options):
        p = GenomeManager()
        p.regenerate_low_tree_red(
            options.split_trees_dir, options.ref_tree,options.red_file)


    def split_reference_tree(self, options):
        p = TreeManager(options.ref_tree, options.tax,
                        options.rank_to_split, options.msa, options.domain)
        make_sure_path_exists(options.out_dir)
        p.split_tree(options.metadata,options.original_log,options.out_dir)

    def merge_logs(self, options):
        merge_logs(options.input_log,options.pruned_tree,options.output_log)

    def unroot_tree(self,options):

        unroot(options.input_tree, options.output_tree)

    def hatchet_wf(self,options):
        make_sure_path_exists(options.out_dir)
        high_level_directory = os.path.join(options.out_dir,'backbone')
        make_sure_path_exists(high_level_directory)
        g = GenomeManager()
        self.logger.info('High level genome picking....')
        g.pick_one_genome(options.ref_tree,options.metadata, options.msa, options.tax,
                          options.domain,options.original_log,
                          options.rank_of_interest,  high_level_directory)


        pruned_tree = os.path.join(high_level_directory, "gtdb_pruned.tree")
        backbone_red_value_file = os.path.join(high_level_directory,'backbone_red_value.tsv')

        self.logger.info('....')
        g.regenerate_red_values(
            options.ref_tree, pruned_tree, options.red_file, backbone_red_value_file)

        class_level_directory = os.path.join(options.out_dir,'class_level')

        t = TreeManager(options.ref_tree, options.tax,
                        options.rank_to_split, options.msa, options.domain)
        make_sure_path_exists(class_level_directory)
        t.split_tree(options.metadata,options.original_log,class_level_directory)

        g.regenerate_low_tree_red(class_level_directory, options.ref_tree, options.red_file)

        to_copy_directory = os.path.join(options.out_dir, 'to_copy')

        package_results(high_level_directory,backbone_red_value_file,class_level_directory,to_copy_directory)

    def parse_options(self, options):
        """Parse user options and call the correct pipeline(s)"""
        if options.subparser_name == 'pick':
            self.pick_genomes(options)
        elif options.subparser_name == 'red':
            self.regenerate_red_values(options)
        elif options.subparser_name == 'split':
            self.split_reference_tree(options)
        elif options.subparser_name == 'red_low':
            self.regenerate_low_tree_red(options)
        elif options.subparser_name == 'merge_logs':
            self.merge_logs(options)
        elif options.subparser_name == 'unroot':
            self.unroot_tree(options)
        elif options.subparser_name == 'hatchet_wf':
            self.hatchet_wf(options)

        else:
            self.logger.error('Unknown command: ' +
                              options.subparser_name + '\n')
            sys.exit()

        return 0
