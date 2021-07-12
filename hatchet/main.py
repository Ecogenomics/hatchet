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


import sys
import logging

from hatchet.biolib_lite.common import make_sure_path_exists
from hatchet.tools import merge_logs

from hatchet.genomemanager import GenomeManager
from hatchet.treemanager import TreeManager


class OptionsParser():
    def __init__(self):
        """Initialization"""

        self.logger = logging.getLogger('timestamp')

    def pick_genomes(self, options):
        make_sure_path_exists(options.out_dir)
        p = GenomeManager()
        p.pick_one_genome(options.ref_tree, options.msa, options.tax,
                          options.domain, options.rank_of_interest, options.out_dir)

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
        p.split_tree(options.out_dir)

    def merge_logs(self, options):
        merge_logs(options.input_log,options.pruned_tree,options.output_log)

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
        else:
            self.logger.error('Unknown command: ' +
                              options.subparser_name + '\n')
            sys.exit()

        return 0
