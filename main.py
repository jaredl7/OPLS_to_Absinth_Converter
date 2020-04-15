#!/usr/bin/env python
import utils
from parameterize import *
from argparse import ArgumentParser


def main():
    parser = ArgumentParser()
    parser.add_argument('-c', '--config', help='The YAML config file to use.', default='config.yaml')
    parser.add_argument('-oa', '--opls-absinth', help='Whether or not to choose the OPLS ABSINTH parameter file for analysis.', action='store_true')
    parser.add_argument('-q', '--quiet', help='Whether or not to output error messages.', action='store_true')
    args = parser.parse_args()

    config = utils.load_config(args.config)

    opls_topology_file          = config['opls-topology-file']
    opls_parameters_file        = config['opls-parameters-file']
    absinth_parameters_file     = config['abs-charmm-pka-file']
    if args.opls_absinth:
        absinth_parameters_file = config['abs-opls-pka-file']
    patched_residue_names_file  = config['patched-residue-names']

    check_cutoff_biotypes_internal(absinth_parameters_file, opls_topology_file, opls_parameters_file, args.quiet)
    add_OPLS_patched_residues(absinth_parameters_file, opls_topology_file, opls_parameters_file, patched_residue_names_file)


if __name__ == '__main__':
    main()
