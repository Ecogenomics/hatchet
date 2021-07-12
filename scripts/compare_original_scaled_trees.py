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

__prog_name__ = 'compare_original_scaled_trees.py'
__prog_desc__ = 'Comapre branch length of the same subtree once it has been rescaled with FasttreeMP.'

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
import dendropy

from collections import defaultdict


class Comparator(object):
    def __init__(self):
        pass

    def run_comparison(self, original_tree, fitted_tree, gtdb_taxonomy,output_file):
        orig_tree = dendropy.Tree.get_from_path(original_tree,
                                           schema='newick',
                                           rooting='force-rooted',
                                           preserve_underscores=True)
        for nd in orig_tree.preorder_node_iter():
            # get left and right taxa that define this node
            taxa = list(nd.preorder_iter(lambda n: n.is_leaf()))
            nd.two_leaves = (taxa[0].taxon.label,taxa[-1].taxon.label)

        gtdb_taxonomy = self.read(gtdb_taxonomy)

        dict_original_length = {}
        for ed in orig_tree.preorder_edge_iter():
            # get left and right taxa that define this node
            child_node = ed.head_node.two_leaves
            if ed.tail_node is not None:
                parent_node = ed.tail_node.two_leaves
                dict_original_length[(parent_node,child_node)] = ed.length

        fit_tree = dendropy.Tree.get_from_path(fitted_tree,
                schema = 'newick',
                rooting = 'force-rooted',
                preserve_underscores = True)
        for nd in fit_tree.preorder_node_iter():
            # get left and right taxa that define this node
            taxa = list(nd.preorder_iter(lambda n: n.is_leaf()))
            nd.two_leaves = (taxa[0].taxon.label,taxa[-1].taxon.label)

        dict_fitted_length = {}
        for ed in fit_tree.preorder_edge_iter():
            # get left and right taxa that define this node
            child_node = ed.head_node.two_leaves
            if ed.tail_node is not None:
                parent_node = ed.tail_node.two_leaves
                dict_fitted_length[(parent_node,child_node)] = ed.length

        count=0
        outfile = open(output_file,'w')
        outfile.write('Edge\tparent_node\tchild_node\thigher_tax\toriginal_length\tnew_length\tlength_diff(abs)\n')
        for k,v in dict_fitted_length.items():
            tax1 = gtdb_taxonomy.get(k[0][0])
            tax2 = gtdb_taxonomy.get(k[0][1])
            combine_tax = [x for idx,x in enumerate(tax1) if tax2[idx]==x]
            list1_as_set = set(tax1)


            outfile.write(f'edge_{count}\t{k[0]}\t{k[1]}\t{";".join(combine_tax)}\t'
                          f'{dict_original_length.get(k)}\t{v}\t'
                          f'{abs(dict_original_length.get(k)-v)}\n')

            count += 1
        outfile.close()

    def read(self, taxonomy_file, use_canonical_gid=False):
        """Read Greengenes-style taxonomy file.
        Expected format is:
            <id>\t<taxonomy string>
        where the taxonomy string has the formats:
            d__; c__; o__; f__; g__; s__
        Parameters
        ----------
        taxonomy_file : str
            Greengenes-style taxonomy file.
        Returns
        -------
        dict : d[unique_id] -> [d__<taxon>, ..., s__<taxon>]
            Taxa indexed by unique ids.
        """

        try:
            d = {}
            for row, line in enumerate(open(taxonomy_file)):
                line_split = line.split('\t')
                unique_id = line_split[0]

                if use_canonical_gid:
                    unique_id = self.canonical_gid(unique_id)

                tax_str = line_split[1].rstrip()
                if tax_str[-1] == ';':
                    # remove trailing semicolons which sometimes
                    # appear in Greengenes-style taxonomy files
                    tax_str = tax_str[0:-1]

                d[unique_id] = [x.strip() for x in tax_str.split(';')]
        except:
            raise

        return d


    def canonical_gid(self,gid):
        """Get canonical form of NCBI genome accession.

        Example:
            G005435135 -> G005435135
            GCF_005435135.1 -> G005435135
            GCF_005435135.1_ASM543513v1_genomic -> G005435135
            RS_GCF_005435135.1 -> G005435135
            GB_GCA_005435135.1 -> G005435135
        """

        if gid.startswith('U'):
            return gid

        gid = gid.replace('RS_', '').replace('GB_', '')
        gid = gid.replace('GCA_', 'G').replace('GCF_', 'G')
        if '.' in gid:
            gid = gid[0:gid.find('.')]

        return gid






if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('original_tree',
                        help='Taxonomy from GTDB-Tk with no split trees.')
    parser.add_argument('fitted_tree',
                        help='Taxonomy from GTDB-Tk with split trees.')
    parser.add_argument('gtdb_taxonomy',
                        help='gtdb taxonomy')
    parser.add_argument('output_file',
                        help='High level taxonomy')
    args = parser.parse_args()

    try:
        gp = Comparator()
        gp.run_comparison(args.original_tree, args.fitted_tree,args.gtdb_taxonomy, args.output_file)
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise
