"""
This file is meant to contain code related to validating and updating the CHARMM parameters in the ABSINTH
parameter using using the OPLS parameter and topology files.
"""
from parsers import OPLS_Topology_Parser, OPLS_Parameters_Parser, Absinth_Parameters_Parser
from collections import OrderedDict
from itertools import chain
from pprint import pprint
import copy
import math


ABSINTH_BIOTYPES_CUTOFF = 1208


def parse_files(config):
    # The absinth parser is unique in that its parameters are built into the object
    absinth_parser = Absinth_Parameters_Parser(config['abs-opls-pka-file'])
    absinth_parser.parse()

    opls_topology_parser = OPLS_Topology_Parser(config['opls-topology-file'])
    opls_topology = opls_topology_parser.parse()

    opls_parameters_parser = OPLS_Parameters_Parser(config['opls-parameters-file'])
    opls_parameters = opls_parameters_parser.parse()

    return absinth_parser, opls_topology, opls_parameters


def locate_value_in_namedtuple(collection_of_tuples, fields, search_value):
    found = list()
    for nt in collection_of_tuples:
        values = [getattr(nt, field) for field in fields]
        if search_value in values:
            found.append(nt)
    return found



def validate_atoms(absinth_parser, opls_masses, opls_residues):
    # valence is the number of unique atoms directly bonded to the atom (in absinth parlance)
    for atom in absinth_parser.atoms:
        possible_matches = list()
        for atom_name in opls_masses:
            opls_atom = opls_masses[atom_name]
            diff = abs(opls_atom.mass - atom.atomic_mass)
            if diff < 1:
                possible_matches.append(opls_atom)
        #print(possible_matches, end='\n\n')
        
        counts = list()
        for match in possible_matches:
            for residue_name in opls_residues:
                residue = opls_residues[residue_name]
                atom_names = [residue.atoms[atom].mass_class for atom in residue.atoms]
                print(match, atom_names)


def validate_biotypes(absinth_parser, opls_topology, opls_parameters):
    # `bond_type` actually acts as a lookup table as it references bonds, angles, torsions, and imptors.
    biotype_connections = OrderedDict()
    for biotype in absinth_parser.biotypes:
    
        # validate the LJ type
        found_LJ = None
        for atom in absinth_parser.atoms:
            if biotype.LJ_type == atom.number:
                found_LJ = copy.copy(atom)
                break

        # validate the charge type
        found_charge = None
        for charge in absinth_parser.charges:
            if biotype.charge_type == charge.number:
                found_charge = copy.copy(charge)
                break

        # validate the bonded types - this appears to be difficult since it requires a lot of lookups
        # first find all the bonds that match the current bond_type in the biotype
        bonds       = locate_value_in_namedtuple(absinth_parser.bonded_type_bonds,      'j,k'.split(','),       biotype.bonded_type)
        angles      = locate_value_in_namedtuple(absinth_parser.bonded_type_angles,     'j,k,l'.split(','),     biotype.bonded_type)
        torsions    = locate_value_in_namedtuple(absinth_parser.bonded_type_torsions,   'j,k,l,m'.split(','),   biotype.bonded_type)
        imptors     = locate_value_in_namedtuple(absinth_parser.bonded_type_imptors,    'j,k,l,m'.split(','),   biotype.bonded_type)

        biotype_connections[biotype] = [bonds, angles, torsions, imptors]
    return biotype_connections


def validate(config):
    absinth_parser, opls_topology, opls_parameters = parse_files(config)
    residues = copy.deepcopy(opls_topology.residues)
    residues.update(opls_topology.patched_residues)
    
    #validate_atoms(absinth_parser, opls_topology.masses, residues)
    #validate_biotypes(absinth_parser, opls_topology, opls_parameters)


def check_cutoff_biotypes_internal(config, suppress_matches=True):
    # This is a self-consistent check to verify that the values match and are being read/consulted properly.
    absinth_parser, opls_topology, opls_parameters = parse_files(config)
    residues = copy.deepcopy(opls_topology.residues)
    residues.update(opls_topology.patched_residues)

    biotypes = [bt for bt in list(absinth_parser.biotypes) if bt.number > ABSINTH_BIOTYPES_CUTOFF]

    # OK. So, here's the problem. In order to check these, I need to know that these atom types exist
    atoms = copy.copy(absinth_parser.atoms)
    charges = copy.copy(absinth_parser.charges)

    # we first need to confirm that the existing types are properly stored.
    for biotype in biotypes:
        LJ_type = [atom for atom in atoms if atom.number == biotype.LJ_type]
        charge_type = [charge for charge in charges if charge.number == biotype.charge_type]
        if len(LJ_type) > 1:
            raise RuntimeError('There are multiple matches for the LJ type in biotype "%d"' % biotype.number)
        if len(charge_type) > 1:
            raise RuntimeError('There are multiple matches for the charge type in biotype "%d"' % biotype.number)

        if not suppress_matches:
            print(biotype)
            print(LJ_type[0])
            print(charge_type[0], end='\n\n')

        if len(LJ_type) == 0:
            print('[ ERROR ] Lennard-Jones (AtomType) not found for biotype %d.' % biotype.number)
            print(biotype, end='\n\n')
            print()

        if len(charge_type) == 0:
            print('[ ERROR ] ChargeType not found for biotype %d.' % biotype.number)
            print(biotype)
            print()