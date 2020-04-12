"""
This file is meant to contain code related to validating and updating the CHARMM parameters in the ABSINTH
parameter using using the OPLS parameter and topology files.
"""
from parsers import OPLS_Topology_Parser, OPLS_Parameters_Parser, Absinth_Parameters_Parser
from pprint import pprint


def parse_files(config):
    opls_top = OPLS_Topology_Parser(config['opls-topology-file'])
    masses, residues, patched_residues = opls_top.parse()

    opls_parser = OPLS_Parameters_Parser(config['opls-parameters-file'])
    opls_params = opls_parser.parse()

    absinth_parser = Absinth_Parameters_Parser(config['abs-charmm-pka-file'])
    absinth_parser.parse()

