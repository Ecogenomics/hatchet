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

__prog_name__ = 'pplacer_lower_rank.py'
__prog_desc__ = 'DESCRIPTION'

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2020'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'p.chaumeil@uq.edu.au'
__status__ = 'Development'

import os
import sys
import argparse
import subprocess
import tempfile
import re
import dendropy

from operator import itemgetter
from external.fastani import FastANI


class PplacerLogger(object):
    """Helper class for writing pplacer output."""

    def __init__(self, fh):
        """Initialise the class.

        Parameters
        ----------
        fh : BinaryIO
            The file to write to .
        """
        self.fh = fh

    def _disp_progress(self, line):
        """Calculates the progress and writes it to stdout.

        Parameters
        ----------
        line : str
            The line passed from pplacer stdout.
        """
        if not line.startswith('working on '):
            sys.stdout.write(
                '\rInitialising pplacer [{}]'.format(line[0:50].center(50)))
            sys.stdout.flush()
        else:
            re_hits = re.search(r'\((\d+)\/(\d+)\)', line)
            current = int(re_hits.group(1))
            total = int(re_hits.group(2))
            sys.stdout.write('\r{}'.format(self._get_progress_str(current,
                                                                  total)))
            sys.stdout.flush()

    def _get_progress_str(self, current, total):
        """Determines the format of the genomes % string.

        Parameters
        ----------
        current : int
            The current number of genomes which have been placed.
        total : int
            The total number of genomes which are to be placed.

        Returns
        -------
        out : str
            A string formatted to show the progress of placement.
        """
        width = 50
        bar = str()
        prop = float(current) / total
        bar += '#' * int(prop * width)
        bar += '-' * (width - len(bar))
        return 'Placing genomes |{}| {}/{} ({:.2f}%)'.format(bar, current,
                                                             total, prop * 100)

    def read(self, line):
        """Reads a line and writes the progress to stdout and the file.

        Parameters
        ----------
        line : str
            A line returned from Prodigal stdout.
        """
        self.fh.write(line)
        line = line.strip()
        self._disp_progress(line)


class Ranker(object):
    def __init__(self, tree_mapping, higher_taxonomy, msa_query, path_to_packages):
        """Initialization."""
        self.tree_map = self._parse_tree_mapping(tree_mapping)
        self.order_map = self._parse_order_mapping(higher_taxonomy)
        self.msa_file = msa_query
        self.af_threshold = 0.65
        self.msa_q = self.read_fasta(msa_query)
        self.gtdb_taxonomy = self.read_taxonomy(
            '/srv/home/uqpchaum/playground/split_gtdbtk_tree/release89/taxonomy/gtdb_taxonomy.tsv')
        self.path_to_packages = path_to_packages
        self.order_rank = ["d__", "p__", "c__", "o__", 'f__', 'g__', 's__']
        self.marker_dict = {"d__": 0.00, "p__": 0.345219220513, "c__": 0.47695132228,
                            "o__": 0.628606202411, "f__": 0.77327890943, "g__": 0.934420004348}

        self.species_radius = self.parse_radius_file()

    def _formatnote(self, sorted_dict, labels):
        """Format the note field by concatenating all information in a sorted dictionary

        Parameters
        ----------
        sorted_dict : sorted dictionary listing reference genomes, ani and alignment fraction for a specific user genome
                    (genomeid, {ani: value, af: value})
        labels : array of label that are removed from the note field

        Returns
        -------
        string
            note field

        """
        note_list = []
        for element in sorted_dict:
            if element[0] not in labels:
                note_str = "{}, {}, {}, {}, {}".format(element[0],
                                                       self.gtdb_taxonomy.get(
                                                           self.add_ncbi_prefix(element[0]))[6],
                                                       self.species_radius.get(
                                                           element[0]),
                                                       round(
                                                           element[1].get('ani'), 2),
                                                       element[1].get('af'))
                note_list.append(note_str)
        return note_list

    def aa_percent_msa(self, aa_string):
        aa_len = sum([1 for c in aa_string if c.isalpha()])
        aa_perc = float(aa_len) / len(aa_string)
        return round(aa_perc * 100, 2)

    def add_ncbi_prefix(self, refname):
        if refname.startswith("GCF_"):
            return "RS_" + refname
        elif refname.startswith("GCA_"):
            return "GB_" + refname
        else:
            return refname

    def read_taxonomy(self, taxonomy_file):
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
            with open(taxonomy_file, 'r') as f:
                for row, line in enumerate(f.readlines()):
                    line_split = line.split('\t')
                    unique_id = line_split[0]

                    tax_str = line_split[1].rstrip()
                    if tax_str[-1] == ';':
                        # remove trailing semicolons which sometimes
                        # appear in Greengenes-style taxonomy files
                        tax_str = tax_str[0:-1]

                    d[unique_id] = [x.strip() for x in tax_str.split(';')]
        except:
            print(
                'Failed to parse taxonomy file on line %d' % (row + 1))
            raise

        return d

    def parse_radius_file(self):
        results = {}
        with open('/srv/home/uqpchaum/playground/split_gtdbtk_tree/release89/radii/gtdb_radii.tsv') as f:
            for line in f:
                infos = line.strip().split('\t')
                gid = infos[1]
                if infos[1].startswith('GB_') or infos[1].startswith('RS_'):
                    gid = gid[3:]
                results[gid] = float(infos[2])
        return results

    def is_float(self, s):
        """Check if a string can be converted to a float.

        Parameters
        ----------
        s : str
            String to evaluate.

        Returns
        -------
        boolean
            True if string can be converted, else False.
        """
        try:
            float(s)
        except ValueError:
            return False

        return True

    def parse_label(self, label):
        """Parse a Newick label which may contain a support value, taxon, and/or auxiliary information.

        Parameters
        ----------
        label : str
            Internal label in a Newick tree.

        Returns
        -------
        float
            Support value specified by label, or None
        str
            Taxon specified by label, or None
        str
            Auxiliary information, on None
        """

        support = None
        taxon = None
        auxiliary_info = None

        if label:
            label = label.strip()
            if '|' in label:
                label, auxiliary_info = label.split('|')

            if ':' in label:
                support, taxon = label.split(':')
                support = float(support)
            else:
                if self.is_float(label):
                    support = float(label)
                elif label != '':
                    taxon = label

        return support, taxon, auxiliary_info

    def make_sure_path_exists(self, path):
        """Create directory if it does not exist.

        Parameters
        ----------
        path : str
            The path to the directory which should be created.

        Returns
        -------
        bool
            True if the path exists.

        Raises
        ------
        BioLibIOException
            If an error was encountered while creating the directory.
        """
        if not path:
            # lack of a path qualifier is acceptable as this
            # simply specifies the current directory
            return True
        elif os.path.isdir(path):
            return True

        try:
            os.makedirs(path)
            return True
        except OSError:
            logger = logging.getLogger('timestamp')
            logger.error('Specified path could not be created: ' + path)
            raise BioLibIOException(
                'Specified path could not be created: ' + path)

    def read_fasta(self, fasta_file, keep_annotation=False):
        """Read sequences from fasta file.

        Parameters
        ----------
        fasta_file : str
            Name of fasta file to read.
        keep_annotation : boolean
            Determine is sequence id should contain annotation.

        Returns
        -------
        dict : dict[seq_id] -> seq
            Sequences indexed by sequence id.
        """

        if not os.path.exists(fasta_file):
            raise InputFileError('Input file %s does not exist.' % fasta_file)

        if os.stat(fasta_file).st_size == 0:
            return {}

        try:

            if fasta_file.endswith('.gz'):
                file_f, file_mode = gzip.open, 'rt'
            else:
                file_f, file_mode = open, 'r'

            seqs = {}
            with file_f(fasta_file, file_mode) as f:

                for line in f.readlines():
                    # skip blank lines
                    if not line.strip():
                        continue

                    if line[0] == '>':
                        if keep_annotation:
                            seq_id = line[1:-1]
                        else:
                            seq_id = line[1:].split(None, 1)[0]

                        seqs[seq_id] = []
                    else:
                        seqs[seq_id].append(line.strip())

            for seq_id, seq in seqs.items():
                seqs[seq_id] = ''.join(seq).replace(' ', '')
        except:
            print(traceback.format_exc())
            print()
            print("[Error] Failed to process sequence file: " + fasta_file)
            sys.exit(1)

        return seqs

    def _parse_tree_mapping(self, tree_map_file):
        result = {}
        with open(tree_map_file) as tmf:
            for line in tmf:
                order, tree_idx = line.strip().split('\t')
                result[order] = tree_idx
        return result

    def standardise_taxonomy(self, taxstring, marker_set=None):
        """Create a 7 rank taxonomy string from an incomplete taxonomy string

        Parameters
        ----------
        taxstring : str
            incomplete taxonomy string
        marker_set : str
            The marker set to use.

        Returns
        -------
        string
            7 rank taxonomy string.
        """
        # return taxstring
        taxlist = taxstring.split(";")
        while '' in taxlist:
            taxlist.remove('')
        #======================================================================
        # if marker_set == 'bac120':
        #     taxlist.insert(0, 'd__Bacteria')
        # if marker_set == 'ar122':
        #     taxlist.insert(0, 'd__Archaea')
        #======================================================================
        taxlist.extend(self.order_rank[len(taxlist):])
        new_taxstring = ";".join(taxlist)
        return new_taxstring

    def _parse_order_mapping(self, taxonomy):
        result = {}
        with open(taxonomy) as tf:
            for line in tf:
                infos = line.split('\t')
                order_name = infos[1].split(';')[3]
                if len(order_name) > 3:
                    result.setdefault(self.tree_map.get(
                        order_name), []).append(infos[0])
        return result

    def tog(self, pplacer_json_out, tree_file):
        """ Convert the pplacer json output into a newick tree.

        Args:
            pplacer_json_out (str): The path to the output of pplacer.
            tree_file (str): The path to output the newick file to.

        Raises:
            TogException: If a non-zero exit code is returned, or the tree file
                          isn't output.
        """

        args = ['guppy', 'tog', '-o', tree_file, pplacer_json_out]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        proc_out, proc_err = proc.communicate()

        if proc.returncode != 0:
            self.logger.error('An error was encountered while running tog.')
            raise TogException(proc_err)

        if not os.path.isfile(pplacer_json_out):
            self.logger.error('tog returned a zero exit code but no output '
                              'file was generated.')
            raise TogException

    def _get_pplacer_taxonomy(self, out_dir, prefix, marker_set_id, user_msa_file, tree):
        """Parse the pplacer tree and write the partial taxonomy for each user genome based on their placements

        Parameters
        ----------
        out_dir : output directory
        prefix : desired prefix for output files
        marker_set_id : bacterial or archaeal id (bac120 or ar122)
        user_msa_file : msa file listing all user genomes for a certain domain
        tree : pplacer tree including the user genomes

        Returns
        -------
        dictionary[genome_label]=pplacer_taxonomy

        """

        out_root = os.path.join(out_dir, 'classify', 'intermediate_results')
        self.make_sure_path_exists(out_root)
        result = {}

        out_pplacer = os.path.join(
            out_dir, '{}.bac120.classification_pplacer.tsv'.format(prefix))

        # We get the pplacer taxonomy for comparison
        with open(out_pplacer, 'w') as pplaceout:
            user_genome_ids = set(self.read_fasta(user_msa_file).keys())
            for leaf in tree.leaf_node_iter():
                if leaf.taxon.label in user_genome_ids:
                    taxa = []

                    cur_node = leaf
                    while cur_node.parent_node:
                        _support, taxon, _aux_info = self.parse_label(
                            cur_node.label)
                        if taxon:
                            for t in taxon.split(';')[::-1]:
                                taxa.append(t.strip())
                        cur_node = cur_node.parent_node
                    taxa.append('d__Bacteria')
                    taxa_str = ';'.join(taxa[::-1])
                    pplaceout.write('{}\t{}\n'.format(
                        leaf.taxon.label, self.standardise_taxonomy(taxa_str, marker_set_id)))
                    result[leaf.taxon.label] = self.standardise_taxonomy(
                        taxa_str, marker_set_id)
        return result

    def _verify_genome_id(self, genome_id):
        """Ensure genome ID will be valid in Newick tree.

        Parameters
        ----------
        genome_id : str
            The string representing the genome identifier.

        Returns
        -------
        bool
            True if the genome identifier is legal.

        Raises
        ------
        GenomeNameInvalid
            If the genome identifier contains illegal characters.
        """

        invalid_chars = set('()[],;=')
        if any((c in invalid_chars) for c in genome_id):
            self.logger.error('Invalid genome ID: %s' % genome_id)
            self.logger.error(
                'The following characters are invalid: %s' % ' '.join(invalid_chars))
            raise GenomeNameInvalid('Invalid genome ID: {}'.format(genome_id))
        return True

    def _get_fastani_genome_path(self, fastani_verification, genomes, tree_file):
        """Generates a queue of comparisons to be made and the paths to
        the corresponding genome id."""
        dict_compare, dict_paths = dict(), dict()

        for qry_node, qry_dict in fastani_verification.items():
            user_label = qry_node.taxon.label
            dict_paths[user_label] = genomes[user_label]
            dict_compare[user_label] = set()
            for node in qry_dict.get('potential_g'):
                leafnode = node[0]
                shortleaf = leafnode.taxon.label
                if leafnode.taxon.label.startswith('GB_') or leafnode.taxon.label.startswith('RS_'):
                    shortleaf = leafnode.taxon.label[3:]
                ref_path = os.path.join(
                    '/srv/home/uqpchaum/playground/split_gtdbtk_tree/release89/fastani/database', shortleaf + "_genomic.fna.gz")

                dict_compare[user_label].add(shortleaf)
                dict_paths[shortleaf] = ref_path

        return dict_compare, dict_paths

    def _genomes_to_process(self, batchfile):
        """Get genomes to process.

        Parameters
        ----------
        genome_dir : str
            Directory containing genomes.
        batchfile : str
            File describing genomes.
        extension : str
            Extension of files to process.

        Returns
        -------
        genomic_files : d[genome_id] -> FASTA file
            Map of genomes to their genomic FASTA files.
        """

        genomic_files = {}
        if batchfile:
            with open(batchfile, "r") as fh:
                for line_no, line in enumerate(fh):
                    line_split = line.strip().split("\t")
                    if line_split[0] == '':
                        continue  # blank line

                    if len(line_split) != 2:
                        self.logger.error(
                            'Batch file must contain exactly 2 columns.')
                        raise GenomeBatchfileMalformed

                    genome_file, genome_id = line_split
                    self._verify_genome_id(genome_id)

                    if genome_file is None or genome_file == '':
                        raise GTDBTkExit(
                            'Missing genome file on line %d.' % (line_no + 1))
                    elif genome_id is None or genome_id == '':
                        raise GTDBTkExit(
                            'Missing genome ID on line %d.' % (line_no + 1))
                    elif genome_id in genomic_files:
                        raise GTDBTkExit(
                            'Genome ID %s appears multiple times.' % genome_id)
                    if genome_file in genomic_files.values():
                        self.logger.warning(
                            'Genome file appears multiple times: %s' % genome_file)

                    genomic_files[genome_id] = genome_file

        # Check that the prefix is valid and the path exists
        invalid_paths = list()
        for genome_key, genome_path in genomic_files.items():
            if genome_key.startswith("RS_") or genome_key.startswith("GB_") \
                    or genome_key.startswith("UBA"):
                self.logger.error("Submitted genomes start with the same prefix"
                                  " (RS_,GB_,UBA) as reference genomes in"
                                  " GTDB-Tk. This will cause issues for"
                                  " downstream analysis.")
                raise GTDBTkExit

            if not os.path.isfile(genome_path):
                invalid_paths.append((genome_key, genome_path))

        # Report on any invalid paths
        if len(invalid_paths) > 0:
            raise GTDBTkExit(
                'There are  paths in the batchfile which do not exist, see gtdb.warnings.log')

        return genomic_files

    def _sort_fastani_results(self, fastani_verification, pplacer_taxonomy_dict,
                              all_fastani_dict, msa_dict, percent_multihit_dict,
                              bac_ar_diff, summaryfout):
        """Format the note field by concatenating all information in a sorted dictionary

        Parameters
        ----------
        fastani_verification : dictionary listing the potential genomes associated with a user genome d[user_genome] = {"potential_g": [
                                    (potential_genome_in_same_genus,patristic distance)], "pplacer_g": genome_of_reference_selected_by_pplacer(if any)}
        all_fastani_dict : dictionary listing the fastani ANI for each user genomes against the potential genomes d[user_genome]={ref_genome1:{"af":af,"ani":ani},ref_genome2:{"af":af,"ani":ani}}
        summaryfout: output file 

        Returns
        -------
        classified_user_genomes: list of genomes where FastANI and Placement in the reference tree have predicted a taxonomy
        unclassified_user_genomes: dictionary of genomes where FastANI and Placement in the reference tree have not  predicted a taxonomy

        """
        classified_user_genomes = []
        unclassified_user_genomes = {}
        for userleaf, potential_nodes in fastani_verification.items():
            summary_list = [None] * 19

            notes = []
            if userleaf.taxon.label in percent_multihit_dict:
                notes.append('Genome has more than {}% of markers with multiple hits'.format(
                    percent_multihit_dict.get(userleaf.taxon.label)))
            if userleaf.taxon.label in bac_ar_diff:
                notes.append('Genome domain questionable ( {}% Bacterial, {}% Archaeal)'.format(
                    bac_ar_diff.get(userleaf.taxon.label).get('bac120'),
                    bac_ar_diff.get(userleaf.taxon.label).get('ar122')))
            if len(notes) > 0:
                summary_list[18] = ';'.join(notes)

            if potential_nodes.get("pplacer_g"):
                pplacer_leafnode = potential_nodes.get("pplacer_g").taxon.label
                if pplacer_leafnode[0:3] in ['RS_', 'GB_']:
                    pplacer_leafnode = pplacer_leafnode[3:]
                if userleaf.taxon.label in all_fastani_dict:
                    # import IPython; IPython.embed()
                    prefilter_reference_dictionary = {k: v for k, v in
                                                      all_fastani_dict.get(userleaf.taxon.label).items() if (
                                                          v.get('ani') >= self.species_radius.get(k) and v.get(
                                                              'af') >= self.af_threshold)}
                    sorted_dict = sorted(iter(all_fastani_dict.get(
                        userleaf.taxon.label).items()), key=lambda _x_y: (_x_y[1]['ani'], _x_y[1]['af']), reverse=True)
                    sorted_prefilter_dict = sorted(iter(prefilter_reference_dictionary.items()),
                                                   key=lambda _x_y1: (_x_y1[1]['ani'], _x_y1[1]['af']), reverse=True)

                    fastani_matching_reference = None
                    if len(sorted_prefilter_dict) > 0:
                        fastani_matching_reference = sorted_prefilter_dict[0][0]
                        current_ani = all_fastani_dict.get(userleaf.taxon.label).get(
                            fastani_matching_reference).get('ani')
                        current_af = all_fastani_dict.get(userleaf.taxon.label).get(
                            fastani_matching_reference).get('af')

                    taxa_str = ";".join(self.gtdb_taxonomy.get(
                        self.add_ncbi_prefix(pplacer_leafnode)))

                    summary_list[0] = userleaf.taxon.label

                    summary_list[11] = pplacer_taxonomy_dict.get(
                        userleaf.taxon.label)
                    summary_list[12] = 'ANI/Placement'
                    summary_list[15] = self.aa_percent_msa(
                        msa_dict.get(summary_list[0]))
                    summary_list[16] = 'XX'

                    if fastani_matching_reference is not None:
                        summary_list[2] = fastani_matching_reference
                        summary_list[3] = str(
                            self.species_radius.get(fastani_matching_reference))
                        summary_list[4] = ";".join(self.gtdb_taxonomy.get(
                            self.add_ncbi_prefix(fastani_matching_reference)))
                        summary_list[5] = round(current_ani, 2)
                        summary_list[6] = current_af
                        if pplacer_leafnode == fastani_matching_reference:
                            if taxa_str.endswith("s__"):
                                taxa_str = taxa_str + pplacer_leafnode
                            summary_list[1] = self.standardise_taxonomy(
                                taxa_str)
                            summary_list[7] = summary_list[2]
                            summary_list[8] = summary_list[4]
                            summary_list[9] = summary_list[5]
                            summary_list[10] = summary_list[6]
                            summary_list[13] = 'topological placement and ANI have congruent species assignments'
                            if len(sorted_dict) > 0:
                                other_ref = '; '.join(self._formatnote(
                                    sorted_dict, [fastani_matching_reference]))
                                if len(other_ref) == 0:
                                    summary_list[14] = None
                                else:
                                    summary_list[14] = other_ref

                        else:
                            taxa_str = ";".join(self.gtdb_taxonomy.get(
                                self.add_ncbi_prefix(fastani_matching_reference)))
                            summary_list[1] = self.standardise_taxonomy(
                                taxa_str)
                            summary_list[7] = pplacer_leafnode
                            summary_list[8] = ";".join(self.gtdb_taxonomy.get(
                                self.add_ncbi_prefix(pplacer_leafnode)))
                            if pplacer_leafnode in all_fastani_dict.get(userleaf.taxon.label):
                                summary_list[9] = round(all_fastani_dict.get(
                                    userleaf.taxon.label).get(pplacer_leafnode).get('ani'), 2)
                                summary_list[10] = all_fastani_dict.get(
                                    userleaf.taxon.label).get(pplacer_leafnode).get('af')
                            summary_list[13] = 'topological placement and ANI have incongruent species assignments'
                            summary_list[12] = 'ANI'

                            if len(sorted_dict) > 0:
                                other_ref = '; '.join(self._formatnote(
                                    sorted_dict, [fastani_matching_reference, pplacer_leafnode]))
                                if len(other_ref) == 0:
                                    summary_list[14] = None
                                else:
                                    summary_list[14] = other_ref

                        summaryfout.write("{}\n".format(
                            '\t'.join(['N/A' if x is None else str(x) for x in summary_list])))
                        classified_user_genomes.append(userleaf.taxon.label)
                    else:
                        summary_list[7] = pplacer_leafnode
                        summary_list[8] = ";".join(self.gtdb_taxonomy.get(
                            self.add_ncbi_prefix(pplacer_leafnode)))
                        if pplacer_leafnode in all_fastani_dict.get(userleaf.taxon.label):
                            summary_list[9] = round(all_fastani_dict.get(
                                userleaf.taxon.label).get(pplacer_leafnode).get('ani'), 2)
                            summary_list[10] = all_fastani_dict.get(
                                userleaf.taxon.label).get(pplacer_leafnode).get('af')

                        if len(sorted_dict) > 0:
                            other_ref = '; '.join(self._formatnote(
                                sorted_dict, [pplacer_leafnode]))
                            if len(other_ref) == 0:
                                summary_list[14] = None
                            else:
                                summary_list[14] = other_ref
                        unclassified_user_genomes[userleaf.taxon.label] = summary_list

            elif userleaf.taxon.label in all_fastani_dict:
                # import IPython; IPython.embed()
                prefilter_reference_dictionary = {k: v for k, v in
                                                  all_fastani_dict.get(userleaf.taxon.label).items() if (
                                                      v.get('ani') >= self.species_radius.get(k) and v.get(
                                                          'af') >= self.af_threshold)}
                sorted_dict = sorted(iter(all_fastani_dict.get(
                    userleaf.taxon.label).items()), key=lambda _x_y2: (_x_y2[1]['ani'], _x_y2[1]['af']), reverse=True)
                sorted_prefilter_dict = sorted(iter(prefilter_reference_dictionary.items()),
                                               key=lambda _x_y3: (_x_y3[1]['ani'], _x_y3[1]['af']), reverse=True)

                summary_list[0] = userleaf.taxon.label
                summary_list[11] = pplacer_taxonomy_dict.get(
                    userleaf.taxon.label)
                summary_list[12] = 'ANI/Placement'
                summary_list[15] = self.aa_percent_msa(
                    msa_dict.get(summary_list[0]))
                summary_list[16] = 'XX'

                if len(sorted_prefilter_dict) > 0:
                    fastani_matching_reference = sorted_prefilter_dict[0][0]

                    taxa_str = ";".join(self.gtdb_taxonomy.get(
                        self.add_ncbi_prefix(fastani_matching_reference))[:-1])

                    summary_list[1] = self.standardise_taxonomy(taxa_str)
                    summary_list[2] = fastani_matching_reference
                    summary_list[3] = str(
                        self.species_radius.get(fastani_matching_reference))
                    summary_list[4] = ";".join(self.gtdb_taxonomy.get(
                        self.add_ncbi_prefix(fastani_matching_reference)))
                    current_ani = all_fastani_dict.get(userleaf.taxon.label).get(
                        fastani_matching_reference).get('ani')
                    summary_list[5] = round(current_ani, 2)
                    current_af = all_fastani_dict.get(userleaf.taxon.label).get(
                        fastani_matching_reference).get('af')
                    summary_list[6] = current_af

                    taxa_str = ";".join(self.gtdb_taxonomy.get(
                        self.add_ncbi_prefix(fastani_matching_reference)))
                    summary_list[1] = self.standardise_taxonomy(
                        taxa_str)

                    summary_list[13] = 'topological placement and ANI have incongruent species assignments'
                    if len(sorted_dict) > 0:
                        other_ref = '; '.join(self._formatnote(
                            sorted_dict, [fastani_matching_reference]))
                        if len(other_ref) == 0:
                            summary_list[14] = None
                        else:
                            summary_list[14] = other_ref

                    summaryfout.write("{}\n".format(
                        '\t'.join(['N/A' if x is None else str(x) for x in summary_list])))

                    classified_user_genomes.append(userleaf.taxon.label)
                else:
                    if len(sorted_dict) > 0:
                        other_ref = '; '.join(self._formatnote(
                            sorted_dict, []))
                        if len(other_ref) == 0:
                            summary_list[14] = None
                        else:
                            summary_list[14] = other_ref
                    unclassified_user_genomes[userleaf.taxon.label] = summary_list
        return classified_user_genomes, unclassified_user_genomes

    def _get_redtax(self, list_subnode, closest_rank):
        """
        Provide a taxonomy string to a user genome based on the reference genomes of the same clade.
        If the clade contains multiple reference genomes we are comparing their taxonomies.
        -If all reference genomes have the same taxonomy up to the 'closest rank' ,
        the taxonomy string including the closest rank is returned.
        -If **NOT** all reference genomes have the same taxonomy up to the 'closest rank',
        the taxonomy string **NOT** including the closest rank is returned.

        Parameters
        ----------
        list_subnode : list of leaf nodes including multiple reference genome.
        closest_rank : last rank of the reference taxonomy

        Returns
        -------
        string
            Taxonomy string.

        """

        subtax, multirefrank = self._parse_subnodes(list_subnode, closest_rank)
        # if all orders in the list are the same, the user genomes gets the
        # same order
        if len(set(multirefrank)) == 1:
                # case d
            subtax.append(multirefrank[0])
        else:
            # otherwise it's stored as undefined
            # case a,b
            subtax.append(closest_rank + "undefined")
        return ';'.join(subtax)

    def _parse_subnodes(self, list_subnode, closest_rank):
        subtax = []
        multirefrank = []
        initial_loop = True
        for item in list_subnode:
            # We get the taxonomy of all reference genomes
            if item.startswith('RS_') or item.startswith('GB_') or item.startswith('UBA'):
                taxonomy_from_file = self.gtdb_taxonomy.get(item)
                # we store the selected rank (i.e. order) for each reference
                # genome
                for rank in taxonomy_from_file:
                    if rank.startswith(closest_rank):
                        multirefrank.append(rank)
                        initial_loop = False
                        break
                    elif initial_loop:
                        # The first iteration is used to stored upper level (
                        # i.e. domain,phylum,class )
                        subtax.append(rank)
        return subtax, multirefrank

    def _assign_mrca_red(self, input_tree, red_file, genomes_to_review):
        """Parse the pplacer tree and write the partial taxonomy for each user genome based on their placements

        Parameters
        ----------
        input_tree : pplacer tree
        marker_set_id : bacterial or archeal id (bac120 or ar122)

        Returns
        -------
        tree: pplacer tree with RED value added to nodes of interest

        """

        glist = []

        print('Calculating RED values based on reference tree.')

        # create map from leave labels to tree nodes
        leaf_node_map = {}
        for leaf in input_tree.leaf_node_iter():
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
                        node = input_tree.mrca(taxa=taxa)
                        node.rel_dist = float(red_value)
                        reference_nodes.add(node)
                elif len(labels) == 1:
                    if labels[0] in leaf_node_map:
                        node = leaf_node_map[labels[0]]
                        node.rel_dist = float(red_value)
                        reference_nodes.add(node)
            seed_node = input_tree.seed_node
            seed_node.rel_dist = 0.00
            reference_nodes.add(seed_node)

        # For all leaf nodes that are not reference genomes
        # We only give RED value to added nodes placed on a reference edge ( between a reference parent and a reference child)
        # The new red value for the pplacer node =
        # RED_parent + (RED_child -RED_parent) * ( (pplacer_disttoroot - parent_disttoroot) / (child_disttoroot - parent_disttoroot) )
        reference_pplacer_node = {}
        for nd in input_tree.leaf_nodes():
            if nd not in reference_nodes:
                nd.rel_dist = 1.0
                pplacer_node = nd
                pplacer_parent_node = pplacer_node.parent_node

                while not bool(set(pplacer_node.leaf_nodes()) & reference_nodes):
                    pplacer_node = pplacer_parent_node
                    pplacer_parent_node = pplacer_node.parent_node

                # perform level-order tree search to find first child
                # node that is part of the reference set
                for child in pplacer_node.levelorder_iter():
                    if child in reference_nodes:
                        child_node = child
                        break

                # find first parent node that is part of the reference set
                while not pplacer_parent_node in reference_nodes:
                    pplacer_parent_node = pplacer_parent_node.parent_node

                if pplacer_parent_node == input_tree.seed_node:
                    genomes_to_review.write(nd.taxon.label + '\n')
                    glist.append(nd.taxon.label)
                    continue

                # we go up the tree until we reach pplacer_parent_node
                current_node = child_node.parent_node
                edge_length = child_node.edge_length
                on_pplacer_branch = False
                pplacer_edge_length = 0

                while current_node != pplacer_parent_node:
                    if on_pplacer_branch or current_node == pplacer_node:
                        on_pplacer_branch = True
                        pplacer_edge_length += current_node.edge_length
                    edge_length += current_node.edge_length
                    current_node = current_node.parent_node

                ratio = pplacer_edge_length / edge_length

                branch_rel_dist = child_node.rel_dist - pplacer_parent_node.rel_dist
                branch_rel_dist = pplacer_parent_node.rel_dist + branch_rel_dist * ratio

                pplacer_node.rel_dist = branch_rel_dist

                #print(pplacer_node, branch_rel_dist, nd.taxon)
        return (input_tree, glist)

    def run(self):
        for k in range(1, 14):
            v = self.order_map.get(str(k))
            genomes = self._genomes_to_process(
                '/srv/home/uqpchaum/playground/split_gtdbtk_tree/clusters_new_genomes/tmp_dir/cluster_{}/genomes.derep.lst'.format(k))

            summaryfout = open(
                '/srv/home/uqpchaum/playground/split_gtdbtk_tree/clusters_new_genomes/tmp_dir/cluster_{}/classifications.tsv'.format(k), 'w')

            genomes_to_review = open(
                '/srv/home/uqpchaum/playground/split_gtdbtk_tree/clusters_new_genomes/tmp_dir/cluster_{}/genomes_to_review.tsv'.format(k), 'w')

            summaryfout.write(
                "user_genome\tclassification\tfastani_reference\tfastani_reference_radius\tfastani_taxonomy\tfastani_ani\tfastani_af\t" +
                "closest_placement_reference\tclosest_placement_taxonomy\tclosest_placement_ani\tclosest_placement_af\tpplacer_taxonomy\t" +
                "classification_method\tnote\tother_related_references(genome_id,species_name,radius,ANI,AF)\taa_percent\ttranslation_table\tred_value\twarnings\n")

            temp_dir = '/srv/home/uqpchaum/playground/split_gtdbtk_tree/clusters_new_genomes/tmp_dir/cluster_{}'.format(
                k)

            self.make_sure_path_exists(temp_dir)
            msa_file = open(os.path.join(temp_dir, 'msa_file.fa'), 'w')
            msa_dict = {}
            for ge in v:
                msa_file.write('>{}\n{}\n'.format(ge, self.msa_q.get(ge)))
                msa_dict[ge] = self.msa_q.get(ge)
            msa_file.close()

            args = ['pplacer', '-m', 'WAG', '-j', str(10), '-c', os.path.join(self.path_to_packages, 'gtdbtk.package.{}.refpkg'.format(k)), '-o',
                    os.path.join(temp_dir, 'pplacer.bac120.json'), os.path.join(temp_dir, 'msa_file.fa')]
            print(' '.join(args))

            #==================================================================
            # proc = subprocess.Popen(
            #     args, stdout=subprocess.PIPE, encoding='utf-8')
            #==================================================================

            #==================================================================
            # proc = subprocess.Popen(
            #     args, stdout=subprocess.PIPE)
            # pplacer_out = os.path.join(temp_dir, 'pplacer.bac120.out')
            # with open(pplacer_out, 'w') as fh:
            #     pplacer_logger = PplacerLogger(fh)
            #     while True:
            #         line = proc.stdout.readline()
            #         if not line:
            #             sys.stdout.write('\n')
            #             break
            #         pplacer_logger.read(line)
            # proc.wait()
            #==================================================================

            self.tog(os.path.join(temp_dir, 'pplacer.bac120.json'),
                     os.path.join(temp_dir, 'pplacer.bac120.tree'))

            # get taxonomic classification of each user genome
            tree_file = os.path.join(temp_dir, 'pplacer.bac120.tree')
            tree = dendropy.Tree.get_from_path(tree_file,
                                               schema='newick',
                                               rooting='force-rooted',
                                               preserve_underscores=True)
            pplacer_taxonomy_dict = self._get_pplacer_taxonomy(
                temp_dir, 'gtdbtk', 'bac120', self.msa_file, tree)

            # Classify step
            # Genomes can be classified by using FastANI or RED values
            # We go through all leaves of the tree. if the leaf is a user
            # genome we take its parent node and look at all the leaves
            # for this node.
            all_fastani_dict = {}
            fastani_verification = {}
            number_comparison = 0
            for userleaf in tree.leaf_node_iter():
                # for each user genome, we select the first parent node with a label.
                # if, while going up the tree, we find a node with only one
                # reference genome, we select this reference genome as
                # leaf_reference.
                if userleaf.taxon.label[0:3] not in ['RS_', 'GB_', 'UBA']:

                    par_node = userleaf.parent_node
                    leaf_ref_genome = None
                    leaf_ref_genomes = [subnd for subnd in par_node.leaf_iter(
                    ) if subnd.taxon.label.replace("'", '')[0:3] in ['RS_', 'GB_', 'UBA']]
                    if len(leaf_ref_genomes) == 1:
                        leaf_ref_genome = leaf_ref_genomes[0]

                    _support, parent_taxon, _aux_info = self.parse_label(
                        par_node.label)
                    # while par_node is not None and parent_taxon is empty,
                    # we go up the tree
                    while par_node is not None and not parent_taxon:
                        par_node = par_node.parent_node
                        if leaf_ref_genome is None:
                            leaf_ref_genomes = [subnd for subnd in par_node.leaf_iter(
                            ) if subnd.taxon.label.replace("'", '')[0:3] in ['RS_', 'GB_', 'UBA']]
                            if len(leaf_ref_genomes) == 1:
                                leaf_ref_genome = leaf_ref_genomes[0]
                        _support, parent_taxon, _aux_info = self.parse_label(
                            par_node.label)

                    # if the parent node is at the genus level
                    parent_rank = parent_taxon.split(";")[-1]
                    if parent_rank.startswith('g__'):
                        # we get all the reference genomes under this genus
                        list_subnode_initials = [subnd.taxon.label.replace(
                            "'", '')[0:3] for subnd in par_node.leaf_iter()]
                        if (list_subnode_initials.count('RS_') + list_subnode_initials.count(
                                'GB_') + list_subnode_initials.count('UBA')) < 1:
                            raise Exception(
                                "There is no reference genomes under '{}'".format('parent_rank'))
                        else:
                            dict_dist_refgenomes = {}
                            list_ref_genomes = [subnd for subnd in par_node.leaf_iter(
                            ) if subnd.taxon.label.replace("'", '')[0:3] in ['RS_', 'GB_', 'UBA']]
                            # we pick the first 100 genomes closest (patristic distance) to the
                            # user genome under the same genus
                            for ref_genome in list_ref_genomes:
                                taxon_labels = [
                                    userleaf.taxon.label, ref_genome.taxon.label]
                                mrca = tree.mrca(taxon_labels=taxon_labels)
                                # the following command is faster than
                                # calculating the patristic distance
                                dict_dist_refgenomes[ref_genome] = (userleaf.distance_from_root(
                                ) - mrca.distance_from_root()) + (
                                    ref_genome.distance_from_root() - mrca.distance_from_root())
                            sorted_l = sorted(
                                iter(dict_dist_refgenomes.items()), key=itemgetter(1))
                            sorted_l = sorted_l[0:100]
                            number_comparison += len(sorted_l)
                            fastani_verification[userleaf] = {
                                "potential_g": sorted_l, "pplacer_g": leaf_ref_genome}
                    else:
                        if leaf_ref_genome:
                            fastani_verification[userleaf] = {"potential_g": [
                                (leaf_ref_genome, 0.0)], "pplacer_g": leaf_ref_genome}

            # we run a fastani comparison for each user genomes against the
            # selected genomes in the same genus
            if len(fastani_verification) > 0:
                fastani = FastANI(cpus=70, force_single=True)

                d_ani_compare, d_paths = self._get_fastani_genome_path(
                    fastani_verification, genomes, tree_file)
                all_fastani_dict = fastani.run(d_ani_compare, d_paths)

            percent_multihit_dict = {}
            bac_ar_diff = {}

            classified_user_genomes, unclassified_user_genomes = self._sort_fastani_results(
                fastani_verification, pplacer_taxonomy_dict, all_fastani_dict, msa_dict, percent_multihit_dict,
                bac_ar_diff, summaryfout)

            print('{0} genome(s) have been classified using FastANI and pplacer.'.format(
                len(classified_user_genomes)))

            # If Fastani can't select a taxonomy for a genome, we use RED
            # distances

            tree_to_process, glist = self._assign_mrca_red(
                tree, '/srv/home/uqpchaum/playground/split_gtdbtk_tree/release89/mrca_red/gtdbtk_r89_bac120.tsv', genomes_to_review)

            user_genome_ids = set(self.read_fasta(self.msa_file).keys())
            # we remove ids already classified with FastANI
            user_genome_ids = user_genome_ids.difference(
                set(classified_user_genomes))
            for leaf in tree_to_process.leaf_node_iter():
                if leaf.taxon.label in glist:
                    continue
                if leaf.taxon.label in user_genome_ids:
                    # In some cases , pplacer can associate 2 user genomes
                    # on the same parent node so we need to go up the tree
                    # to find a node with a reference genome as leaf.
                    cur_node = leaf.parent_node
                    list_subnode_initials = [subnd.taxon.label.replace(
                        "'", '')[0:3] for subnd in cur_node.leaf_iter()]
                    while 'RS_' not in list_subnode_initials and 'GB_' not in list_subnode_initials and 'UBA' not in list_subnode_initials:
                        cur_node = cur_node.parent_node
                        list_subnode_initials = [subnd.taxon.label.replace(
                            "'", '')[0:3] for subnd in cur_node.leaf_iter()]

                    current_rel_list = cur_node.rel_dist

                    parent_taxon_node = cur_node.parent_node
                    _support, parent_taxon, _aux_info = self.parse_label(
                        parent_taxon_node.label)

                    while parent_taxon_node is not None and not parent_taxon:
                        parent_taxon_node = parent_taxon_node.parent_node
                        _support, parent_taxon, _aux_info = self.parse_label(
                            parent_taxon_node.label)

                    # is the node represent multiple ranks, we select the lowest one
                    # i.e. if node is p__A;c__B;o__C we pick o__
                    parent_rank = parent_taxon.split(";")[-1][0:3]
                    if not hasattr(parent_taxon_node, 'rel_dist'):
                        genomes_to_review.write(
                            '{}\n'.format(leaf.taxon.label))
                        continue
                    parent_rel_dist = parent_taxon_node.rel_dist

                    debug_info = [leaf.taxon.label, parent_rank,
                                  parent_rel_dist, '', '', '', '']

                    child_taxons = []
                    closest_rank = None
                    detection = "taxonomic novelty determined using RED"
                    # if the genome is not placed between the genus and
                    # specie ranks
                    if parent_rank != 'g__':
                        # we select the child rank (if parent_rank = 'c__'
                        # child rank will be 'o__)'
                        child_rk = self.order_rank[self.order_rank.index(
                            parent_rank) + 1]

                        # get all reference genomes under the current node
                        list_subnode = [childnd.taxon.label.replace("'", '') for childnd in cur_node.leaf_iter(
                        ) if childnd.taxon.label[0:3] in ['RS_', 'UBA', 'GB_']]

                        # get all names for the child rank
                        list_ranks = [self.gtdb_taxonomy.get(
                            name)[self.order_rank.index(child_rk)] for name in list_subnode]

                        # if there is just one rank name
                        if len(set(list_ranks)) == 1:
                            for subranknd in cur_node.preorder_iter():
                                _support, subranknd_taxon, _aux_info = self.parse_label(
                                    subranknd.label)
                                if subranknd.is_internal() and subranknd_taxon is not None and subranknd_taxon.startswith(
                                        child_rk):
                                    child_taxons = subranknd_taxon.split(
                                        ";")
                                    child_taxon_node = subranknd
                                    child_rel_dist = child_taxon_node.rel_dist
                                    break
                        else:
                            # case 2a and 2b
                            closest_rank = parent_rank
                            detection = "taxonomic classification fully defined by topology"
                    else:
                        # case 1a
                        closest_rank = parent_rank
                        detection = "taxonomic classification fully defined by topology"

                    # case 1b
                    if len(child_taxons) == 0 and closest_rank is None:
                        list_leaves = [childnd.taxon.label.replace("'", '') for childnd in cur_node.leaf_iter(
                        ) if childnd.taxon.label[0:3] in ['RS_', 'UBA', 'GB_']]
                        if len(list_leaves) != 1:
                            list_subrank = []
                            for leaf_subrank in list_leaves:
                                list_subrank.append(self.gtdb_taxonomy.get(
                                    leaf_subrank)[self.order_rank.index(parent_rank) + 1])
                            if len(set(list_subrank)) == 1:

                                print(list_leaves)
                                print(list_subrank)
                                print(set(list_subrank))
                                print(leaf.taxon.label)
                                raise Exception(
                                    'There should be only one leaf.')
                            else:
                                closest_rank = parent_rank
                                detection = "taxonomic classification fully defined by topology"
                        list_leaf_ranks = self.gtdb_taxonomy.get(
                            list_leaves[0])[self.order_rank.index(child_rk):-1]  # We remove the species name
                        for leaf_taxon in reversed(list_leaf_ranks):
                            if leaf_taxon == list_leaf_ranks[0]:
                                if abs(current_rel_list - self.marker_dict.get(leaf_taxon[:3])) < abs(
                                        current_rel_list - self.marker_dict.get(parent_rank)):
                                    closest_rank = leaf_taxon[:3]
                                    debug_info[3] = leaf_taxon
                                    debug_info[5] = 'case 1b - III'
                                    break
                            else:
                                pchildrank = list_leaf_ranks[list_leaf_ranks.index(
                                    leaf_taxon) - 1]
                                if abs(current_rel_list - self.marker_dict.get(leaf_taxon[:3])) < abs(
                                        current_rel_list - self.marker_dict.get(pchildrank[:3])):
                                    closest_rank = leaf_taxon[:3]
                                    debug_info[1] = pchildrank
                                    debug_info[2] = 1.0
                                    debug_info[3] = leaf_taxon
                                    debug_info[5] = 'case 1b - II'
                                    break
                        if closest_rank is None:
                            closest_rank = parent_rank
                            debug_info[3] = list_leaf_ranks[0]
                            debug_info[5] = 'case 1b - IV'

                    # if there is multiple ranks on the child node (i.e genome between p__Nitrospirae and c__Nitrospiria;o__Nitrospirales;f__Nitropiraceae)
                    # we loop through the list of rank from f_ to c_ rank
                    for child_taxon in reversed(child_taxons):
                        # if lower rank is c__Nitropiria
                        if child_taxon == child_taxons[0]:
                            if (abs(current_rel_list - self.marker_dict.get(child_taxon[:3])) < abs(
                                    child_rel_dist - self.marker_dict.get(child_taxon[:3])) and
                                    abs(current_rel_list - self.marker_dict.get(child_taxon[:3])) < abs(
                                        current_rel_list - self.marker_dict.get(parent_rank))):
                                debug_info[3] = ';'.join(child_taxons)
                                debug_info[4] = child_rel_dist
                                debug_info[5] = 'case 3b - II'
                                closest_rank = child_taxon[:3]
                            elif closest_rank is None:
                                closest_rank = parent_rank
                                debug_info[3] = ';'.join(child_taxons)
                                debug_info[4] = child_rel_dist
                                debug_info[5] = 'case 3b - III'
                        else:
                            pchildrank = child_taxons[child_taxons.index(
                                child_taxon) - 1]
                            if (abs(current_rel_list - self.marker_dict.get(child_taxon[:3])) < abs(
                                    current_rel_list - self.marker_dict.get(pchildrank[:3])) and
                                    abs(current_rel_list - self.marker_dict.get(child_taxon[:3])) < abs(
                                        child_rel_dist - self.marker_dict.get(child_taxon[:3]))):
                                closest_rank = child_taxon
                                debug_info[3] = ';'.join(child_taxons)
                                debug_info[4] = child_rel_dist
                                debug_info[5] = 'case 3b - I'
                                break

                    # case 1b
                    if closest_rank is None:
                        raise Exception('closest rank is None')

                    debug_info[6] = closest_rank

                    list_subnode = [subnd.taxon.label.replace(
                        "'", '') for subnd in cur_node.leaf_iter()]
                    red_taxonomy = self._get_redtax(
                        list_subnode, closest_rank)

                    del debug_info[0]

                    summary_list = [None] * 19
                    if leaf.taxon.label in unclassified_user_genomes:
                        summary_list = unclassified_user_genomes.get(
                            leaf.taxon.label)
                        if summary_list[13] == '':
                            summary_list[13] = None
                    summary_list[0] = leaf.taxon.label
                    summary_list[1] = self.standardise_taxonomy(
                        red_taxonomy)
                    summary_list[11] = pplacer_taxonomy_dict.get(
                        leaf.taxon.label)
                    summary_list[12] = 'Placement'
                    summary_list[13] = detection
                    summary_list[15] = self.aa_percent_msa(
                        msa_dict.get(summary_list[0]))
                    summary_list[16] = 'XX'
                    summary_list[17] = current_rel_list

                    notes = []
                    if summary_list[0] in percent_multihit_dict:
                        notes.append('Genome has more than {}% of markers with multiple hits'.format(
                            percent_multihit_dict.get(summary_list[0])))
                    if summary_list[0] in bac_ar_diff:
                        notes.append('Genome domain questionable ( {}% Bacterial, {}% Archaeal)'.format(
                            bac_ar_diff.get(summary_list[0]).get('bac120'),
                            bac_ar_diff.get(summary_list[0]).get('ar122')))

                    if len(notes) > 0:
                        summary_list[18] = ';'.join(notes)
                    summaryfout.write("{0}\n".format(
                        '\t'.join(['N/A' if x is None else str(x) for x in summary_list])))

            summaryfout.close()
            genomes_to_review.close()


if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--tree_mapping', help='')
    parser.add_argument('--higher_taxonomy', help='')
    parser.add_argument('--msa_query', help='')
    parser.add_argument('--path_to_packages')
    args = parser.parse_args()

    try:
        ranker = Ranker(args.tree_mapping,
                        args.higher_taxonomy, args.msa_query, args.path_to_packages)
        ranker.run()
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise
