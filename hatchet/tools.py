
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

