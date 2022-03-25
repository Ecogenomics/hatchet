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

__prog_name__ = 'copy_and_rename_bac_summary.py'
__prog_desc__ = 'find every gtdb bac120 summary file and move them to another directory, while changing their name'

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2021'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'p.chaumeil@uq.edu.au'
__status__ = 'Development'

import argparse
import sys


sys.setrecursionlimit(15000)

class Mover(object):
    def __init__(self):
        self.rank_order = ['d', 'p', 'c', 'o', 'f', 'g', 's']

    def run(self):



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
        tt = Mover()
        tt.run()
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise




