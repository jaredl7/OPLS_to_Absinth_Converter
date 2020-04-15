"""
This file is meant to contain code related to validating and updating the CHARMM parameters in the ABSINTH
parameter using using the OPLS parameter and topology files.
"""
from parsers import OPLS_Topology_Parser, OPLS_Parameters_Parser, Absinth_Parameters_Parser
from common_objects import BioType, AtomType
from collections import OrderedDict
from itertools import chain
from pprint import pprint
from utils import load_yaml_file
import copy
import math
import re


ABSINTH_BIOTYPES_CUTOFF = 1208


def parse_files(absinth_parameter_file, opls_topology_file, opls_parameter_file):
    # The absinth parser is unique in that its parameters are built into the object
    absinth_parser = Absinth_Parameters_Parser(absinth_parameter_file)
    absinth_parser.parse()

    opls_topology_parser = OPLS_Topology_Parser(opls_topology_file)
    opls_topology = opls_topology_parser.parse()

    opls_parameters_parser = OPLS_Parameters_Parser(opls_parameter_file)
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
    
    validate_atoms(absinth_parser, opls_topology.masses, residues)
    validate_biotypes(absinth_parser, opls_topology, opls_parameters)


def check_cutoff_biotypes_internal(absinth_parameters_file, opls_topology_file, opls_parameters_file, suppress_matches=True):
    # This is a self-consistent check to verify that the values match and are being read/consulted properly.
    absinth_parser, opls_topology, opls_parameters = parse_files(absinth_parameters_file, opls_topology_file, opls_parameters_file)
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
        elif len(LJ_type) == 0:
            print('[ ERROR ] Lennard-Jones (AtomType) not found for biotype %d.' % biotype.number)
            print(biotype, end='\n\n')
            print()

        if len(charge_type) > 1:
            raise RuntimeError('There are multiple matches for the charge type in biotype "%d"' % biotype.number)
        elif len(charge_type) == 0:
            print('[ ERROR ] ChargeType not found for biotype %d.' % biotype.number)
            print(biotype)
            print()

        if not suppress_matches:
            print(biotype)
            if len(LJ_type) == 0:
                print(LJ_type)
            else:
                print(LJ_type[0])

            if len(charge_type) == 0:
                print(charge_type, end='\n\n')
            else:
                print(charge_type[0], end='\n\n')


def search_biotypes(patched_residues, absinth_parser):
    available_LJ_types = list(set(atom.name for atom in absinth_parser.atoms))
    regex = re.compile('([A-Za-z]+)([0-9]+)')
    numbers = '0123456789'
    for residue_name in patched_residues:
        residue = patched_residues[residue_name]
        
        # For this to work we need to 
        for atom in residue.atoms:
            matches = [match for match in regex.split(atom) if len(match) > 0]
            if len(matches) > 1:
                if matches[-1] in numbers:
                    print(atom)
            #print(atom, matches)


def determine_atomic_number(search_atom):
    atomic_numbers = {'MW': 0, 'H': 1, 'C': 6, 'N': 7, 'O': 8, 'NA': 11, 'P': 15, 'S': 16, 'CL': 17, 'K': 19, 'BR': 35, 'I': 53, 'CS': 55}
    regex_primary = re.compile('([A-Za-z]+)([0-9]+)')
    regex_alternative = re.compile('([0-9]+)([A-Za-z]+)')

    regex = regex_primary
    if search_atom[0].isnumeric():
        regex = regex_alternative
    split_search_string = [match for match in regex.split(search_atom) if len(match) > 0]
    
    search = None
    result = None
    # this is for the case where the search_atom is something like "C" or "N"
    if len(split_search_string) == 1:
        search = split_search_string[0]  # this is guaranteed to be a non-number

    # This is for the case where the search atom is of the format "NXX" or "XXN"
    elif len(split_search_string) == 2:
        search = split_search_string[0]
        if search.isnumeric():
            search = split_search_string[-1]

    query = search[::]
    if query in atomic_numbers:
        result = atomic_numbers[query]
    else:
        query = search[0]
        if query in atomic_numbers:
            result = atomic_numbers[query]
    return result


def determine_atom_valence(residue, search_atom):
    # To figure out the number of atoms, we have to check everything: bonds, angles, and dihedrals
    valence = 0
    connections = list()
    for bond in residue.bonds:
        if search_atom in bond:
            connections.append(tuple(sorted(bond)))

    for angle in residue.angles:
        if search_atom in angle:
            index = angle.index(search_atom)
            if index == 0:
                connections.append(tuple(sorted(angle[:2])))
            elif index == 1:
                connections.append(tuple(sorted(angle[:2])))
                connections.append(tuple(sorted(angle[1:])))
            else:
                connections.append(tuple(sorted(angle[1:])))

    for dihedral in residue.dihedrals:
        if search_atom in dihedral:
            index = dihedral.index(search_atom)
            if index == 0:
                connections.append(tuple(sorted(angle[:2])))
            elif index == 1:
                connections.append(tuple(sorted(angle[:2])))
                connections.append(tuple(sorted(angle[1:3])))
            else:
                connections.append(tuple(sorted(angle[1:3])))
                connections.append(tuple(sorted(angle[2:])))
    set_connections = list(set(connections))
    valence = len(set_connections)
    return valence


def add_OPLS_LJ_types(residue, opls_topology, current_atom_type):
    new_atom_types = list()
    for atom_name in residue.atoms:
        possible_LJ_type = atom_name
        trailing_description = ""

        # this is for the case where the first letter is a number
        if possible_LJ_type[0].isnumeric():
            possible_LJ_type = atom_name[1:]

        # now define the attributes for the `LJ_type`
        number = copy.copy(current_atom_type)
        element = copy.copy(possible_LJ_type)
        description = "\"OPLS %s (%s)\"" % (residue.name, atom_name)
        atomic_number = determine_atomic_number(atom_name)

        atom_mass_class = residue.atoms[atom_name].mass_class
        atom_mass = opls_topology.masses[atom_mass_class].mass
        valence = determine_atom_valence(residue, atom_name)

        new_atom_type = AtomType(current_atom_type, atom_name, description, atomic_number, atom_mass, valence)
        new_atom_types.append(new_atom_type)
        current_atom_type += 1
        # print(new_atom_type)
    return new_atom_types, current_atom_type
        


def add_OPLS_patched_residues(absinth_parameters_file, opls_topology_file, opls_parameters_file, patched_residue_names_file):
    # This is the directive: take the patched OPLS residues, and add them to the biotypes.
    absinth_parser, opls_topology, opls_parameters = parse_files(absinth_parameters_file, opls_topology_file, opls_parameters_file)
    patched_names = load_yaml_file(patched_residue_names_file)['friendly-names']
    current_atom_type = absinth_parser.atoms[-1].number + 1

    # OK - so as I understand it, this has to work in the following manner:
    # To create a biotype, we need to define:
    # 1) A LJ-type
    # 2) A charge type
    # 3) A bond type
    #
    # The bond type is mainly a lookup table. While the others are definitions for how things should be put together.
    all_new_atom_types = copy.copy(absinth_parser.atoms)
    for residue_name in opls_topology.patched_residues:
        residue = opls_topology.patched_residues[residue_name]

        # the easiest way into this is to create new biotypes, and ignore everything that
        # exists previously. I'll write another version that accounts for the previous.
        new_atom_types, current_atom_type = add_OPLS_LJ_types(residue, opls_topology, current_atom_type)
        all_new_atom_types += new_atom_types  # update the atom_types / LJ_types


