import subprocess

import dendropy

from hatchet.biolib_lite.common import check_file_exists
from hatchet.biolib_lite.newick import parse_label


def purge_reload(command, shell_command_file, cmd):
    shell_command_file.write('module purge;')
    shell_command_file.write('module load {};'.format(command))
    shell_command_file.write('{};'.format(cmd))


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

def remove_character(decorated_tree,char):
    # first get all lines from file
    with open(decorated_tree, 'r') as f:
        lines = f.readlines()
    # remove character
    lines = [line.replace(char, '') for line in lines]
    # finally, write lines in the file
    with open(decorated_tree, 'w') as f:
        f.writelines(lines)

def execute(cmd):
    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)

def unroot(rooted_tree, unrooted_tree):
    # get taxonomic classification of each user genome
    tree = dendropy.Tree.get_from_path(rooted_tree,
                                       schema='newick',
                                       rooting='force-rooted',
                                       preserve_underscores=True)
    tree.deroot()
    tree.write_to_path(unrooted_tree,
                       schema='newick',
                       suppress_rooting=True,
                       unquoted_underscores=True)
    return True

def prune(selflogger,input_tree, taxa_to_retain_file, output_tree):
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

    selflogger.info('Reading input tree.')
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
    selflogger.info('Identifying taxa to retain.')
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

    selflogger.info('Identified %d extant taxa to retain in tree.' % len(taxa_in_tree))

    if taxa_to_retain:
        selflogger.warning('Failed to identify %d taxa: %s' % (len(taxa_to_retain), ','.join(taxa_to_retain)))

    # prune tree
    selflogger.info('Pruning tree.')
    tree.retain_taxa(taxa_in_tree)

    # write out results
    selflogger.info('Writing output tree.')
    tree.write_to_path(output_tree,
                       schema='newick',
                       suppress_rooting=True,
                       unquoted_underscores=True)

