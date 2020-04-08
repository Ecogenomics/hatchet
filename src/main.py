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
from collections import defaultdict

from biolib.common import check_file_exists, make_sure_path_exists, check_dir_exists

from src.genomemanager import GenomeManager
from src.treemanager import TreeManager


class OptionsParser():
    def __init__(self):
        """Initialization"""

        self.logger = logging.getLogger('timestamp')

    def pick_genomes(self, options):
        make_sure_path_exists(options.output_dir)
        p = GenomeManager()
        p.pick_one_genome(options.tree, options.msa, options.taxonomy,
                          options.domain, options.rank_of_interest, options.output_dir)

    def regenerate_red_values(self, options):
        p = GenomeManager()
        p.regenerate_red_values(
            options.raw_tree, options.pruned_tree, options.red_file, options.output)

    def regenerate_low_tree_red(self, options):
        p = GenomeManager()
        p.regenerate_low_tree_red(
            options.split_trees_dir, options.reference_tree,options.red_file)

    def split_reference_tree(self, options):
        p = TreeManager(options.reference_tree, options.gtdb_taxonomy,
                        options.rank_to_split, options.msa_reference, options.domain)
        p.split_tree(options.output_dir)

    def parse_options(self, options):
        """Parse user options and call the correct pipeline(s)"""
        if options.subparser_name == 'pick':
            self.pick_genomes(options)
        elif options.subparser_name == 'recreate_red':
            self.regenerate_red_values(options)
        elif options.subparser_name == 'split':
            self.split_reference_tree(options)
        elif options.subparser_name == 'red_low':
            self.regenerate_low_tree_red(options)
        else:
            self.logger.error('Unknown command: ' +
                              options.subparser_name + '\n')
            sys.exit()

        return 0
