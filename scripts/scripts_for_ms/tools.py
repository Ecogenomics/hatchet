import errno
import gzip
import logging
import os
import sys
import traceback

import dendropy


def make_sure_path_exists(path):
    """Create directory if it does not exist."""

    if not path:
        # lack of a path qualifier is acceptable as this
        # simply specifies the current directory
        return

    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            logger = logging.getLogger('timestamp')
            logger.error('Specified path could not be created: ' + path + '\n')
            sys.exit()

def parse_label(label):
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
            label, auxiliary_info = label.split('|', 1)

        if ':' in label:
            support, taxon = label.split(':')
            support = float(support)
        else:
            if is_float(label):
                support = float(label)
            elif label != '':
                taxon = label

    return support, taxon, auxiliary_info

def is_float(s):
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

def read_fasta(fasta_file, keep_annotation=False):
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
        raise OSError('Input file %s does not exist.' % fasta_file)

    if os.stat(fasta_file).st_size == 0:
        return {}

    try:
        open_file = open
        if fasta_file.endswith('.gz'):
            open_file = gzip.open

        seqs = {}
        for line in open_file(fasta_file, 'rt'):
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



def check_file_exists(input_file):
    """Check if file exists."""
    if not os.path.exists(input_file) or not os.path.isfile(input_file):
        logger = logging.getLogger('timestamp')
        logger.error('Input file does not exists: ' + input_file + '\n')
        sys.exit()

def prune(input_tree, taxa_to_retain_file, output_tree):
    """Prune tree.

    Parameters
    ----------
    input_tree : str
      File containing Newick tree to rerooted.
    taxa_to_retain_file : str
      File specifying taxa to retain.
    output_tree : str
      Name of file for rerooted tree.
    """

    check_file_exists(input_tree)
    check_file_exists(taxa_to_retain_file)

    logging.info('Reading input tree.')
    tree = dendropy.Tree.get_from_path(input_tree,
                                       schema='newick',
                                       rooting='force-rooted',
                                       preserve_underscores=True)

    # read taxa to retain
    taxa_to_retain = set()
    for line in open(taxa_to_retain_file):
        if line[0] == '#' or not line.strip():
            continue

        line_split = line.strip().split('\t')
        taxa_to_retain.add(line_split[0])

    # find taxa to retain
    logging.info('Identifying taxa to retain.')
    taxa_in_tree = set()
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            support, taxon, _auxiliary_info = parse_label(node.taxon.label)
            if taxon in taxa_to_retain:
                taxa_in_tree.add(node.taxon)
                taxa_to_retain.remove(taxon)
        else:
            support, taxon, _auxiliary_info = parse_label(node.label)
            if taxon in taxa_to_retain:
                for leaf in node.leaf_iter():
                    taxa_in_tree.add(leaf.taxon)
                taxa_to_retain.remove(taxon)

        # check if all outgroup taxa have been identified
        if not taxa_to_retain:
            break

    logging.info('Identified %d extant taxa to retain in tree.' % len(taxa_in_tree))

    if taxa_to_retain:
        logging.warning('Failed to identify %d taxa: %s' % (len(taxa_to_retain), ','.join(taxa_to_retain)))

    # prune tree
    logging.info('Pruning tree.')
    tree.retain_taxa(taxa_in_tree)

    # write out results
    logging.info('Writing output tree.')
    tree.write_to_path(output_tree,
                       schema='newick',
                       suppress_rooting=True,
                       unquoted_underscores=True)


def regenerate_red_values(raw_tree, pruned_trees, red_file, output):
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

    write_rd(pruned_tree, output, unpruned_tree)

def write_rd(pruned_tree, output_rd_file, unpruned_tree):
    """Write out relative divergences for each node."""

    fout = open(output_rd_file, 'w')
    for n in pruned_tree.preorder_node_iter():
        if n.is_leaf():
            fout.write('%s\t%f\n' % (n.taxon.label, 1.00))
            #fout.write('%s\t%f\n' % (n.taxon.label, n.rel_dist))
        else:
            # get left and right taxa that define this node
            taxa = list(n.preorder_iter(lambda n: n.is_leaf()))
            # get rel_dist of this node in the original tree
            reldist_node = unpruned_tree.mrca(
                taxon_labels=[taxa[0].taxon.label, taxa[-1].taxon.label])
            fout.write('%s|%s\t%f\n' %
                       (taxa[0].taxon.label, taxa[-1].taxon.label, reldist_node.rel_dist))

    fout.close()


def merge_logs(log_fitting_in,non_fitted_tree,log_fitting_out):
    print(log_fitting_in,non_fitted_tree,log_fitting_out)
    output_file = open(log_fitting_out, 'w')
    with open(non_fitted_tree, 'r') as fistripped:
        stripped_tree = fistripped.read().replace('\n', '')
    last_item = 0
    with open(log_fitting_in) as logf:
        for line in logf:
            if line.startswith('ML_Lengths'):
                last_item += 1

    nber_iter = 0
    with open(log_fitting_in) as logf:
        for line in logf:
            if line.startswith('ML_Lengths'):
                nber_iter += 1
                if nber_iter == last_item:
                    output_file.write(f"ML_Lengths\t{stripped_tree}\n")
                else:
                    output_file.write(line)
            else:
                output_file.write(line)