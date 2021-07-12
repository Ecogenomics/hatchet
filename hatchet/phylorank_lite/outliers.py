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

import logging
import os
import random
import sys
import copy
from collections import defaultdict, namedtuple

import dendropy
from biolib.plots.abstract_plot import AbstractPlot
from biolib.taxonomy import Taxonomy
from numpy import (median as np_median,
                   array as np_array,
                   arange as np_arange,
                   percentile as np_percentile,
                   ones_like as np_ones_like,
                   histogram as np_histogram)

from phylorank.common import (read_taxa_file,
                              filter_taxa_for_dist_inference,
                              get_phyla_lineages)
from phylorank.newick import parse_label
from rel_dist import RelativeDistance
from phylorank.viral_taxonomy import (VIRAL_RANK_LABELS,
                                      translate_viral_taxonomy,
                                      translate_viral_tree,
                                      rev_translate_output_file,
                                      read_viral_taxonomy_from_tree)
import mpld3


class AxisReplacer(mpld3.plugins.PluginBase):
    """Replaces the y axis labels with the provided list"""

    JAVASCRIPT = """
    mpld3.register_plugin("axisreplacer", AxisReplacer);
    AxisReplacer.prototype = Object.create(mpld3.Plugin.prototype);
    AxisReplacer.prototype.constructor = AxisReplacer;
    AxisReplacer.prototype.requiredProps = ["data"];
    AxisReplacer.prototype.defaultProps = {}
    function AxisReplacer(fig, props){
        mpld3.Plugin.call(this, fig, props);
    };

    AxisReplacer.prototype.draw = function(){
        let data = this.props.data;
        let parent = document.getElementsByClassName("mpld3-yaxis")[0];
        let gTicks = parent.getElementsByTagName("g");
        for (let i=0; i < gTicks.length; i++) {
            let curTick = gTicks[i];
            let curText = curTick.getElementsByTagName("text")[0];
            curText.innerHTML = data[i];
        }        
    };
    """

    def __init__(self, linedata):
        self.dict_ = {"type": "axisreplacer",
                      "data": linedata}


class Outliers(AbstractPlot):
    """Identify outliers based on relative distances.

    A named group is considered an outlier if it falls
    far from the median of all other groups at the same
    rank. Since the relative distance of a group from
    the root depending on the rooting, the input tree
    is rooted on all phyla. Results are reported for
    each of these rooting and a consensus file produced
    indicating the median distance over all rootings.
    """

    def __init__(self, skip_mpld3=False, dpi=96, output_dir=None):
        """Initialize."""

        self.logger = logging.getLogger()

        self.skip_mpld3 = skip_mpld3
        if not self.skip_mpld3:
            import mpld3

        Options = namedtuple('Options', 'width height tick_font_size label_font_size dpi')
        options = Options(5, 4, 12, 12, 96)

        AbstractPlot.__init__(self, options)

        self.poly_color = (89.0 / 255, 89.0 / 255, 89.0 / 255)
        self.near_mono_color = (255.0 / 255, 188.0 / 255, 121.0 / 255)
        self.mono_color = (95.0 / 255, 158.0 / 255, 209.0 / 255)

        self.median_color = (0.0 / 255, 107.0 / 255, 164.0 / 255)

        self.dpi = dpi
        self.output_dir = output_dir

    def root_with_outgroup(self, input_tree, taxonomy, outgroup_taxa):
        """Reroot the tree using the given outgroup.

        Parameters
        ----------
        input_tree : Dendropy Tree
          Tree to rerooted.
        taxonomy : dict
          Taxonomy for taxa.
        outgroup_taxa : str
          Desired outgroup.

        Returns
        -------
        Dendropy Tree
            Deep-copy of original tree rerooted on outgroup.
        """

        new_tree = input_tree.clone()

        outgroup = set()
        for genome_id, taxa in taxonomy.items():
            if outgroup_taxa in taxa:
                outgroup.add(genome_id)
        self.logger.info('Identifying %d genomes in the outgroup.' % len(outgroup))

        outgroup_in_tree = set()
        ingroup_in_tree = set()
        for n in new_tree.leaf_node_iter():
            if n.taxon.label in outgroup:
                outgroup_in_tree.add(n.taxon)
            else:
                ingroup_in_tree.add(n)
        self.logger.info('Identified %d outgroup taxa in the tree.' % len(outgroup_in_tree))

        if len(outgroup_in_tree) == 0:
            self.logger.warning('No outgroup taxa identified in the tree.')
            self.logger.warning('Tree was not rerooted.')
            sys.exit(0)

        # There is a complication here. We wish to find the MRCA of the outgroup
        # taxa. Finding the MRCA requires a rooted tree and we have no guarantee
        # that the tree isn't currently rooted within the outgroup clade. There is
        # also no way to identify a node that is guaranteed to be outside the outgroup
        # clade. As such, the tree is randomly rooted on a leaf node not in the outgroup.
        # This random rerooting is performed until the MRCA does not spans all taxa in
        # the tree.

        leaves_in_tree = sum([1 for _ in new_tree.leaf_node_iter()])
        while True:
            rnd_ingroup_leaf = random.sample(ingroup_in_tree, 1)[0]
            new_tree.reroot_at_edge(rnd_ingroup_leaf.edge,
                                    length1=0.5 * rnd_ingroup_leaf.edge_length,
                                    length2=0.5 * rnd_ingroup_leaf.edge_length)

            mrca = new_tree.mrca(taxa=outgroup_in_tree)
            leaves_in_mrca = sum([1 for _ in mrca.leaf_iter()])
            if leaves_in_mrca != leaves_in_tree:
                break

        if leaves_in_mrca != len(outgroup_in_tree):
            self.logger.info('Outgroup is not monophyletic. Tree will be rerooted at the MRCA of the outgroup.')
            self.logger.info('The outgroup consisted of %d taxa, while the MRCA has %d leaf nodes.' % (
                len(outgroup_in_tree), leaves_in_mrca))
            if leaves_in_mrca == leaves_in_tree:
                self.logger.warning('The MRCA spans all taxa in the tree.')
                self.logger.warning('This indicating the selected outgroup is likely polyphyletic in the current tree.')
                self.logger.warning('Polyphyletic outgroups are not suitable for rooting. Try another outgroup.')
        else:
            self.logger.info('Outgroup is monophyletic.')

        if mrca.edge_length is None:
            self.logger.info('Tree appears to already be rooted on this outgroup.')
        else:
            self.logger.info('Rerooting tree.')
            new_tree.reroot_at_edge(mrca.edge,
                                    length1=0.5 * mrca.edge_length,
                                    length2=0.5 * mrca.edge_length)

        return new_tree

    def _median_outlier_file(self,
                             rel_dists,
                             taxa_for_dist_inference,
                             gtdb_parent_ranks,
                             viral,
                             output_file,
                             median_rank_file):
        """Identify outliers relative to the median of rank distributions.

        Parameters
        ----------
        rel_dists: d[rank_index][taxon] -> relative divergence
            Relative divergence of taxa at each rank.
        taxa_for_dist_inference : iterable
            Taxa to considered when inferring distributions.
        gtdb_parent_ranks: d[taxon] -> string indicating parent taxa
            Parent taxa for each taxon.
        output_file : str
            Desired name of output table.
        """

        # determine median relative distance for each rank
        median_rel_dist = {}
        for rank, d in rel_dists.items():
            v = [dist for taxa, dist in d.items() if taxa in taxa_for_dist_inference]
            if len(v) == 0:
                continue

            median_rel_dist[rank] = np_median(v)

        fout_rank = open(median_rank_file, 'w')
        median_str = []
        for rank in sorted(median_rel_dist.keys()):
            if viral:
                rank_label = VIRAL_RANK_LABELS[rank]
            else:
                rank_label = Taxonomy.rank_labels[rank]
            median_str.append('"' + rank_label + '":' + str(median_rel_dist[rank]))
        fout_rank.write('{' + ','.join(median_str) + '}\n')
        fout_rank.close()

        fout = open(output_file, 'w')
        fout.write('Taxa\tGTDB taxonomy\tMedian distance\tMean difference\tClosest rank\tClassification\n')

        for i, rank in enumerate(sorted(rel_dists.keys())):
            for clade_label, dist in rel_dists[rank].items():
                if rank in median_rel_dist:
                    delta = dist - median_rel_dist[rank]
                    closest_rank_dist = 1e10
                    for test_rank, test_median in median_rel_dist.items():
                        abs_dist = abs(dist - test_median)
                        if abs_dist < closest_rank_dist:
                            closest_rank_dist = abs_dist
                            closest_rank = Taxonomy.rank_labels[test_rank]

                    classification = "OK"
                    if delta < -0.2:
                        classification = "very overclassified"
                    elif delta < -0.1:
                        classification = "overclassified"
                    elif delta > 0.2:
                        classification = "very underclassified"
                    elif delta > 0.1:
                        classification = "underclassified"

                    fout.write('%s\t%s\t%.3f\t%.3f\t%s\t%s\n' % (clade_label,
                                                                 ';'.join(gtdb_parent_ranks[clade_label]),
                                                                 dist,
                                                                 delta,
                                                                 closest_rank,
                                                                 classification))
                else:
                    fout.write('%s\t%s\t%.3f\t%.3f\t%s\t%s\n' % (clade_label,
                                                                 ';'.join(gtdb_parent_ranks[clade_label]),
                                                                 dist,
                                                                 -1,
                                                                 'NA',
                                                                 'Insufficent data to calcualte median for rank.'))
        fout.close()

    def taxa_median_rd(self, phylum_rel_dists):
        """Calculate the median relative divergence for each taxon.

        Parameters
        ----------
        phylum_rel_dists: phylum_rel_dists[phylum][rank_index][taxon] -> relative divergences
            Relative divergence of taxon at each rank for different phylum-level rootings.
        """

        medians_for_taxa = defaultdict(lambda: defaultdict(list))
        for p in phylum_rel_dists:
            for rank, d in phylum_rel_dists[p].items():
                for taxon, dist in d.items():
                    medians_for_taxa[rank][taxon].append(dist)

        return medians_for_taxa

    def rank_median_rd(self, phylum_rel_dists, taxa_for_dist_inference):
        """Calculate median relative divergence for each rank.

        Parameters
        ----------
        phylum_rel_dists: phylum_rel_dists[phylum][rank_index][taxon] -> relative divergences
            Relative divergence of taxon at each rank for different phylum-level rootings.
        taxa_for_dist_inference : iterable
            Taxa to considered when inferring distributions.
        """

        medians_for_taxa = self.taxa_median_rd(phylum_rel_dists)

        median_for_rank = {}
        for i, rank in enumerate(sorted(medians_for_taxa.keys())):
            v = [np_median(dists) for taxon, dists in medians_for_taxa[rank].items() if
                 taxon in taxa_for_dist_inference]

            if v:
                median_for_rank[rank] = np_median(v)

        return median_for_rank

    def _median_summary_outlier_file(self, phylum_rel_dists,
                                     taxa_for_dist_inference,
                                     gtdb_parent_ranks,
                                     outlier_table,
                                     rank_file,
                                     verbose_table):
        """Identify outliers relative to the median of rank distributions.

        Parameters
        ----------
        phylum_rel_dists: phylum_rel_dists[phylum][rank_index][taxon] -> relative divergences
            Relative divergence of taxon at each rank for different phylum-level rootings.
        taxa_for_dist_inference : iterable
            Taxa to considered when inferring distributions.
        gtdb_parent_ranks: d[taxon] -> string indicating parent taxa
            Parent taxa for each taxon.
        outlier_table : str
            Desired name of output table.
        rank_file : str
            Desired name of file indicating median relative distance of each rank.
        verbose_table : boolean
            Print additional columns in output table.
        """

        # determine median relative distance for each taxa
        medians_for_taxa = self.taxa_median_rd(phylum_rel_dists)

        # determine median relative distance for each rank
        median_for_rank = self.rank_median_rd(phylum_rel_dists, taxa_for_dist_inference)

        fout_rank = open(rank_file, 'w')
        median_str = []
        for rank in sorted(median_for_rank.keys()):
            median_str.append('"' + Taxonomy.rank_labels[rank] + '":' + str(median_for_rank[rank]))
        fout_rank.write('{' + ','.join(median_str) + '}\n')
        fout_rank.close()

        fout = open(outlier_table, 'w')
        if verbose_table:
            fout.write('Taxa\tGTDB taxonomy\tMedian distance')
            fout.write('\tMedian of rank\tMedian difference')
            fout.write('\tClosest rank\tClassifciation\n')
        else:
            fout.write('Taxa\tGTDB taxonomy\tMedian distance\tMedian difference\tClosest rank\tClassification\n')

        for rank in sorted(median_for_rank.keys()):
            for clade_label, dists in medians_for_taxa[rank].items():
                dists = np_array(dists)

                taxon_median = np_median(dists)
                delta = taxon_median - median_for_rank[rank]

                closest_rank_dist = 1e10
                for test_rank, test_median in median_for_rank.items():
                    abs_dist = abs(taxon_median - test_median)
                    if abs_dist < closest_rank_dist:
                        closest_rank_dist = abs_dist
                        closest_rank = Taxonomy.rank_labels[test_rank]

                classification = "OK"
                if delta < -0.2:
                    classification = "very overclassified"
                elif delta < -0.1:
                    classification = "overclassified"
                elif delta > 0.2:
                    classification = "very underclassified"
                elif delta > 0.1:
                    classification = "underclassified"

                if verbose_table:
                    fout.write('%s\t%s\t%.2f\t%.3f\t%.3f\t%s\t%s\n' % (clade_label,
                                                                       ';'.join(gtdb_parent_ranks[clade_label]),
                                                                       taxon_median,
                                                                       median_for_rank[rank],
                                                                       delta,
                                                                       closest_rank,
                                                                       classification))
                else:
                    fout.write('%s\t%s\t%.3f\t%.3f\t%s\t%s\n' % (clade_label,
                                                                 ';'.join(gtdb_parent_ranks[clade_label]),
                                                                 taxon_median,
                                                                 delta,
                                                                 closest_rank,
                                                                 classification))
        fout.close()

    def rd_fixed_root(self, tree):
        """Scale tree and calculate relative divergence over a single fixed root.

        Parameters
        ----------
        tree : Tree
          Dendropy tree.
        taxa_for_dist_inference : set
          Taxa to use for inference relative divergence distributions.
        """

        # calculate relative distance to taxa
        rd = RelativeDistance()
        rel_dists = rd.rel_dist_to_named_clades(tree)

        # create scaled tree
        rd.decorate_rel_dist(tree)
        for n in tree.preorder_node_iter(lambda n: n != tree.seed_node):
            rd_to_parent = n.rel_dist - n.parent_node.rel_dist
            n.edge_length = rd_to_parent

        return rel_dists


    def median_rd_over_phyla(self,
                             tree,
                             taxa_for_dist_inference,
                             taxonomy=None):
        """Calculate the median relative divergence over all phyla rootings.

        Parameters
        ----------
        tree : Tree
          Dendropy tree.
        taxa_for_dist_inference : set
          Taxa to use for inference relative divergence distributions.
        """

        # read taxonomy from tree
        if not taxonomy:
            taxonomy = Taxonomy().read_from_tree(tree,
                                                 warnings=False)

        # get list of phyla level lineages
        all_phyla = get_phyla_lineages(tree)
        self.logger.info('Identified %d phyla.' % len(all_phyla))

        phyla = [p for p in all_phyla if p in taxa_for_dist_inference]
        self.logger.info('Using %d phyla as rootings for inferring distributions.' % len(phyla))
        if len(phyla) < 2:
            self.logger.error('Rescaling requires at least 2 valid phyla.')
            sys.exit(1)

        # give each node a unique id
        for i, n in enumerate(tree.preorder_node_iter()):
            n.id = i

        # calculate relative divergence for tree rooted on each phylum
        phylum_rel_dists = {}
        rel_node_dists = defaultdict(list)
        rd = RelativeDistance()
        for p in phyla:
            phylum = p.replace('p__', '').replace(' ', '_').lower()
            self.logger.info('Calculating information with rooting on %s.' % phylum.capitalize())

            cur_tree = self.root_with_outgroup(tree, taxonomy, p)

            # calculate relative distance to taxa
            rel_dists = rd.rel_dist_to_named_clades(cur_tree)
            rel_dists.pop(0, None)  # remove results for Domain

            # remove named groups in outgroup
            children = Taxonomy().children(p, taxonomy)
            for r in rel_dists.keys():
                rel_dists[r].pop(p, None)

            for t in children:
                for r in rel_dists.keys():
                    rel_dists[r].pop(t, None)

            phylum_rel_dists[phylum] = rel_dists

            # calculate relative distance to all nodes
            rd.decorate_rel_dist(cur_tree)

            # determine which lineages represents the 'ingroup'
            ingroup_subtree = None
            for c in cur_tree.seed_node.child_node_iter():
                _support, taxon_name, _auxiliary_info = parse_label(c.label)
                if not taxon_name or p not in taxon_name:
                    ingroup_subtree = c
                    break

            # do a preorder traversal of 'ingroup' and record relative divergence to nodes
            for n in ingroup_subtree.preorder_iter():
                rel_node_dists[n.id].append(n.rel_dist)

        return phylum_rel_dists, rel_node_dists

    def _write_rd(self, tree, output_rd_file):
        """Write out relative divergences for each node."""

        fout = open(output_rd_file, 'w')
        for n in tree.preorder_node_iter():
            if n.is_leaf():
                fout.write('%s\t%f\n' % (n.taxon.label, n.rel_dist))
            else:
                # get left and right taxa that define this node
                taxa = list(n.preorder_iter(lambda n: n.is_leaf()))
                fout.write('%s|%s\t%f\n' % (taxa[0].taxon.label, taxa[-1].taxon.label, n.rel_dist))

        fout.close()

    def _write_rd_tree(self, tree, rel_node_dists, output_tree):
        """Write out tree with RED specified at each internal node."""

        # copy tree so node labels aren't changed in original tree
        red_tree = copy.deepcopy(tree)

        for node_id, n in enumerate(red_tree.preorder_node_iter()):
            if n == red_tree.seed_node:
                red = 0
            else:
                red = np_median(rel_node_dists[node_id])

            red_str = "|RED={:.3f}".format(red)
            if n.is_leaf():
                n.taxon.label += red_str
            else:
                if n.label:
                    n.label += red_str
                else:
                    n.label = red_str

        red_tree.write_to_path(output_tree,
                               schema='newick',
                               suppress_rooting=True,
                               unquoted_underscores=True)

    def read_fmeasure(self, fmeasure_table):
        """Read table with F-measure for taxa."""

        fmeasure = {}
        with open(fmeasure_table) as f:
            f.readline()

            for line in f:
                line_split = line.strip().split('\t')

                taxa = line_split[0]
                fmeasure[taxa] = float(line_split[2])

        return fmeasure

    def run(self,
            input_tree,
            taxonomy_file,
            plot_taxa_file,
            plot_dist_taxa_only,
            plot_domain,
            highlight_polyphyly,
            highlight_taxa_file,
            trusted_taxa_file,
            fixed_root,
            min_children,
            min_support,
            mblet,
            fmeasure_table,
            min_fmeasure,
            fmeasure_mono,
            verbose_table):
        """Determine distribution of taxa at each taxonomic rank.

        Parameters
        ----------
        input_tree : str
          Name of input tree.
        taxonomy_file : str
          File with taxonomy strings for each taxa.
        plot_taxa_file : str
          File specifying taxa to plot. Set to None to consider all taxa.
        plot_dist_taxa_only : boolean
          Only plot the taxa used to infer distribution.
        plot_domain : boolean
          Plot domain rank.
        trusted_taxa_file : str
          File specifying trusted taxa to consider when inferring distribution. Set to None to consider all taxa.
        fixed_root : boolean
          Usa a single fixed root to infer outliers.
        min_children : int
          Only consider taxa with at least the specified number of children taxa when inferring distribution.
        min_support : float
          Only consider taxa with at least this level of support when inferring distribution.
        verbose_table : boolean
          Print additional columns in output table.
        """

        # read tree
        self.logger.info('Reading tree.')
        tree = dendropy.Tree.get_from_path(input_tree,
                                           schema='newick',
                                           rooting='force-rooted',
                                           preserve_underscores=True)

        input_tree_name = os.path.splitext(os.path.basename(input_tree))[0]

        # pull taxonomy from tree and file
        self.logger.info('Reading taxonomy.')
        taxonomy = Taxonomy().read(taxonomy_file)

        self.logger.info('Reading taxonomy from tree.')
        tree_taxonomy = Taxonomy().read_from_tree(input_tree,
                                                  warnings=False)

        gtdb_parent_ranks = Taxonomy().parents(tree_taxonomy)


        # determine taxa to be used for inferring distribution
        taxa_for_dist_inference = filter_taxa_for_dist_inference(tree,
                                                                 taxonomy,
                                                                 min_children,
                                                                 min_support)
        self.logger.info('Identified {:,} taxa for use in inferring RED distribution.'.format(
            len(taxa_for_dist_inference)))

        # limit plotted taxa
        if plot_dist_taxa_only:
            taxa_to_plot = taxa_for_dist_inference
        elif plot_taxa_file:
            taxa_to_plot = read_taxa_file(plot_taxa_file)
        else:
            # plot every taxon defined in tree
            taxa_to_plot = set()
            for node in tree.preorder_node_iter():
                support, taxon, _auxiliary_info = parse_label(node.label)
                if taxon:
                    taxon = taxon.split(';')[-1].strip()  # get most specific taxon from compound names
                    # (e.g. p__Armatimonadetes; c__Chthonomonadetes)
                    taxa_to_plot.add(taxon)

            if False:
                # HACK FOR NCBI: only plot taxa with >= 2 taxa
                taxa_to_plot = set()
                for node in tree.preorder_node_iter():
                    if not node.label or node.is_leaf():
                        continue

                    support, taxon, _auxiliary_info = parse_label(node.label)
                    if not taxon:
                        continue
                    taxon = taxon.split(';')[-1].strip()  # get most specific taxon from compound names
                    # (e.g. p__Armatimonadetes; c__Chthonomonadetes)

                    # count number of subordinate children
                    rank_prefix = taxon[0:3]
                    if min_children > 0 and rank_prefix != 's__':
                        child_rank_index = Taxonomy().rank_index[rank_prefix] + 1
                        child_rank_prefix = Taxonomy.rank_prefixes[child_rank_index]
                        subordinate_taxa = set()
                        for leaf in node.leaf_iter():
                            taxa = taxonomy.get(leaf.taxon.label, Taxonomy.rank_prefixes)
                            if len(taxa) > child_rank_index:
                                sub_taxon = taxa[child_rank_index]
                                if sub_taxon != Taxonomy.rank_prefixes[child_rank_index] and sub_taxon.startswith(
                                        child_rank_prefix):
                                    subordinate_taxa.add(sub_taxon)

                        if len(subordinate_taxa) < min_children:
                            continue

                    taxa_to_plot.add(taxon)

        # highlight taxa
        highlight_taxa = set()
        if highlight_taxa_file:
            for line in open(highlight_taxa_file):
                highlight_taxa.add(line.strip().split('\t')[0])

        # check if a single fixed root should be used
        summary_median_outlier_table = os.path.join(self.output_dir, '{}.tsv'.format(input_tree_name))
        distribution_table = os.path.join(self.output_dir, '{}.rank_distribution.tsv'.format(
            input_tree_name))
        summary_median_rank_file = os.path.join(self.output_dir, '{}.dict'.format(input_tree_name))
        phyla_file = os.path.join(self.output_dir, '{}.phyla.tsv'.format(input_tree_name))


        rel_dists = self.rd_fixed_root(tree, taxa_for_dist_inference)

        # create fixed rooting style tables and plots
        self._median_outlier_file(rel_dists,
                                  taxa_for_dist_inference,
                                  gtdb_parent_ranks,
                                  viral,
                                  summary_median_outlier_table,
                                  summary_median_rank_file)

        output_rd_file = os.path.join(self.output_dir, '{}.node_rd.tsv'.format(input_tree_name))
        self._write_rd(tree, output_rd_file)

        output_tree = os.path.join(self.output_dir, '{}.scaled.tree'.format(input_tree_name))
        tree.write_to_path(output_tree,
                           schema='newick',
                           suppress_rooting=True,
                           unquoted_underscores=True)
