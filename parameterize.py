"""
This file is meant to contain code related to validating and updating the CHARMM parameters in the ABSINTH
parameter using using the OPLS parameter and topology files.
"""
from parsers import OPLS_Topology_Parser, OPLS_Parameters_Parser, Absinth_Parameters_Parser
from common_objects import BioType, AtomType, ChargeType, InteractType, ContactType
from common_objects import BondLength, BondAngle, DihedralAngle
from common_objects import BondedTypeBond, BondedTypeAngle, BondedTypeTorsion, BondedTypeImptors
from common_objects import Bond, Angle, Dihedral, Improper
from collections import OrderedDict
from itertools import chain
from pprint import pprint
from utils import load_yaml_file
from tabulate import tabulate
import copy
import math
import re


ABSINTH_BIOTYPES_CUTOFF = 1208
SCALING_CONSTANT = 2.0**(1.0/6.0)


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
    residues = opls_topology.residues[::]
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


def determine_search_string(search_atom):
    regex_primary = re.compile('([A-Za-z]+)([0-9]+)')
    regex_alternative = re.compile('([0-9]+)([A-Za-z]+)')

    regex = regex_primary
    if search_atom[0].isnumeric():
        regex = regex_alternative
    split_search_string = [match for match in regex.split(search_atom) if len(match) > 0]
    
    search = None
    # this is for the case where the search_atom is something like "C" or "N"
    if len(split_search_string) == 1:
        search = split_search_string[0]  # this is guaranteed to be a non-number

    # This is for the case where the search atom is of the format "NXX" or "XXN"
    elif len(split_search_string) == 2:
        search = split_search_string[0]
        if search.isnumeric():
            search = split_search_string[-1]
    return search


def determine_atomic_number(search_atom):
    atomic_numbers = {'MW': 0, 'H': 1, 'C': 6, 'N': 7, 'O': 8, 'NA': 11, 'P': 15, 'S': 16, 'CL': 17, 'K': 19, 'BR': 35, 'I': 53, 'CS': 55}
    search = determine_search_string(search_atom)
    result = None

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
        #print(new_atom_type)
        #print(atom_name, atom_mass_class)
    return new_atom_types, current_atom_type


def add_OPLS_contact_and_interact_types(residue, opls_parameters):
    # TODO: REVISE
    pass


def add_OPLS_bonds(opls_parameters, starting_bond_number, potential_type=1):
    # This one is a little tricky. For the bondlengths (i.e. "bond" in Absinth parlance),
    # the two atoms are not required when defining the bond (that's done in bonded_type_bond).
    # So, all we have to do here is create a new "bond" from the last bond number.

    '''
    # IMPORTANT: If I understand this correctly, the bonds for ABSINTH are completely different than the bonds
    # described in OPLS. OPLS defines two atoms and how they should be paired, while ABSINTH is completely
    # abstract. The best way I can think of this is that `bonded_type_bond` refers to two different
    # bonds which will be built from chemistry and geometry. Enumerating the possible pairs allows for
    # the bonds to be created.
    bond_combinations = OrderedDict()
    for opls_bond in opls_parameters.bonds:
        a = opls_bond.atom_a
        b = opls_bond.atom_b
        pair = (a, b)
        bond = BondLength(starting_bond_number, potential_type, opls_bond.spring_constant, opls_bond.reference_bond_length, 0)  # there is no `c`
        bond_combinations[pair] = bond
        starting_bond_number += 1

        bond = BondLength(starting_bond_number, potential_type, opls_bond.spring_constant, opls_bond.reference_bond_length, 0)  # there is no `c`
        bond_combinations[pair[::-1]] = bond
        starting_bond_number += 1
    return list(bond_combinations.values())    
    '''
    bonds = list()
    for opls_bond in opls_parameters.bonds:
        #print(opls_bond)
        bond = BondLength(starting_bond_number, potential_type, opls_bond.spring_constant, opls_bond.reference_bond_length, 0)  # there is no `c`
        starting_bond_number += 1
        bonds.append(bond)
    return bonds
    

def add_OPLS_angles(opls_parameters, starting_angle_number, potential_type=1):
    angles = list()
    for opls_angle in opls_parameters.angles:
        angle = BondAngle(starting_angle_number, potential_type, opls_angle.theta_force_constant, opls_angle.reference_angle, 0, 0)  # No Urey-Bradley terms
        starting_angle_number += 1
        angles.append(angle)
    return angles


def add_OPLS_dihedrals(opls_parameters, starting_dihedral_number, potential_type=1):
    dihedrals = list()
    for opls_dihedral in opls_parameters.dihedrals:
        dihedral = DihedralAngle(starting_dihedral_number, potential_type, 0, 0, 0, 0, 0, 0, 0)  # No Urey-Bradley terms
        starting_dihedral_number += 1
        dihedrals.append(dihedral)
    return dihedrals


def add_OPLS_bonded_types(opls_topology, opls_parameters, starting_bond_number):
    bonded_type_bonds = list()
    for bt in opls_parameters.bonds:
        #print(bt)
        # figure out what's the name
        starting_bond_number += 1
    return bonded_type_bonds


def add_OPLS_biotypes(residue, opls_topology, current_biotype):
    new_biotypes = list()
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

        current_biotype += 1
    return new_biotypes, current_biotype


def create_atom_mass_class_lookup_table(residues):
    atoms_lookup_table = OrderedDict()
    for residue_name in residues:
        residue = residues[residue_name]
        for atom_name in residue.atoms:
            atom = residue.atoms[atom_name]
            if atom_name not in atoms_lookup_table:
                atoms_lookup_table[atom_name] = [atom.mass_class]
            else:
                atoms_lookup_table[atom_name].append(atom.mass_class)
                atoms_lookup_table[atom_name] = sorted(list(set(atoms_lookup_table[atom_name])))
    return atoms_lookup_table


def create_bond_lookup_table(bonds):
    bond_lookup_table = OrderedDict()
    for bond in bonds:
        key = (bond.atom_a, bond.atom_b)
        bond_lookup_table[key] = bond
    return bond_lookup_table


def find_line_contents(data, search_string):
    line_numbers = list()
    for number, line in enumerate(data):
        if line.startswith(search_string):
            line_numbers.append(number)
    return line_numbers


def output_object_to_list(generic_object, name):
    # note: this is only for namedtuples!
    output = [name]
    for field in generic_object._fields:
        value = getattr(generic_object, field)
        output.append(value)
    return output


def output_objects_to_string(object_container, name):
    objects = list()
    for generic_object in object_container:
        object_as_list = output_object_to_list(generic_object, name)
        objects.append(object_as_list)
    new_lines = tabulate(objects, object_container[0]._fields, tablefmt='plain')
    return new_lines.split('\n')[1:]


def add_OPLS_patched_residues(absinth_parameters_file, opls_topology_file, opls_parameters_file, patched_residue_names_file, savename):
    absinth_parser, opls_topology, opls_parameters = parse_files(absinth_parameters_file, opls_topology_file, opls_parameters_file)
    
    # first read in the parameters file
    with open(absinth_parameters_file) as f:
        absinth_raw = f.read().split('\n')

    # identify the locations of each section - these will determine where the new parameters will be added
    absinth_atom_locations                  = find_line_contents(absinth_raw, 'atom ')
    absinth_contact_locations               = find_line_contents(absinth_raw, 'contact ')
    absinth_interact_locations              = find_line_contents(absinth_raw, 'interact ')
    absinth_radius_locations                = find_line_contents(absinth_raw, 'radius ')
    absinth_charge_locations                = find_line_contents(absinth_raw, 'charge ')
    absinth_fos_locations                   = find_line_contents(absinth_raw, 'fos ')
    absinth_bond_locations                  = find_line_contents(absinth_raw, 'bond ')
    absinth_angle_locations                 = find_line_contents(absinth_raw, 'angle ')
    absinth_torsion_locations               = find_line_contents(absinth_raw, 'torsion ')
    absinth_biotype_locations               = find_line_contents(absinth_raw, 'biotype ')
    absinth_bondedtype_bond_locations       = find_line_contents(absinth_raw, 'bonded_type_bond ')
    absinth_bondedtype_angle_locations      = find_line_contents(absinth_raw, 'bonded_type_angle ')
    absinth_bondedtype_imptors_locations    = find_line_contents(absinth_raw, 'bonded_type_imptors ')
    absinth_bondedtype_torsions_locations   = find_line_contents(absinth_raw, 'bonded_type_torsions ')

    current_atom_type               = absinth_parser.atoms[-1].number + 1
    current_contact_type            = absinth_parser.contacts[-1].i + 1
    current_interact_type           = absinth_parser.interacts[-1].i + 1
    current_charge_type             = absinth_parser.charges[-1].number + 1
    current_biotype                 = absinth_parser.biotypes[-1].number + 1
    current_bonded_type             = list(set(sorted([bt.bond_potential for bt in absinth_parser.bonded_type_bonds])))[-1] + 1
    current_bond_number             = absinth_parser.bonds[-1].identifier + 1
    current_angle_number            = absinth_parser.angles[-1].identifier + 1
    current_dihedral_number         = absinth_parser.torsions[-1].identifier + 1
    current_bondedtype_bond_number  = copy.copy(current_bond_number) + 1

    new_bonds       = add_OPLS_bonds(opls_parameters, current_bond_number)
    new_angles      = add_OPLS_angles(opls_parameters, current_angle_number)
    new_dihedrals   = add_OPLS_dihedrals(opls_parameters, current_dihedral_number)

    nonbonded_atoms = [nb.atom for nb in opls_parameters.nonbonded]
    new_atom_types              = list()
    new_contacts_by_mass_class  = OrderedDict()
    new_contacts_by_number      = OrderedDict()
    new_interacts_by_mass_class = OrderedDict()
    new_interacts_by_number     = OrderedDict()
    new_contacts                = list()
    new_interacts               = list()
    new_charges                 = list()
    new_biotypes                = list()
    # new_bonded_type_bonds       = list()
    # new_bonded_type_angles      = list()
    # new_bonded_type_torsions    = list()
    # new_bonded_type_imptors     = list()
    #
    # The bond type is mainly a lookup table. While the others are definitions for how things should be put together.
    all_new_atom_types          = list()
    all_residue_mass_classes    = list()

    # populate the biotypes, charge types, and interact types
    for residue_name in opls_topology.patched_residues:
        residue = opls_topology.patched_residues[residue_name]
        residue_atom_types, current_atom_type = add_OPLS_LJ_types(residue, opls_topology, current_atom_type)
        new_atom_types += residue_atom_types  # update the atom_types / LJ_types

        # need to add the contact and interact parameters for the atom / LJ_types
        residues_biotypes = list()
        for idx, atom_name in enumerate(residue.atoms):
            atom = residue.atoms[atom_name]
            index = nonbonded_atoms.index(atom.mass_class)

            if atom.mass_class not in new_contacts_by_mass_class:
                new_contacts_by_mass_class[atom.mass_class] = (atom_name, copy.copy(current_contact_type))
                new_contacts_by_number[current_contact_type] = (atom_name, atom.mass_class)
                current_contact_type += 1

                match = new_contacts_by_mass_class[atom.mass_class]
                i = copy.copy(match[-1])
                j = copy.copy(match[-1])
                sigma = 2 * opls_parameters.nonbonded[index].Rmin2 / SCALING_CONSTANT  # Rmin2 means (Rmin / 2)
                contact_type = ContactType(i, j, sigma)
                new_contacts.append(contact_type)

            if atom.mass_class not in new_interacts_by_mass_class:
                new_interacts_by_mass_class[atom.mass_class] = (atom_name, copy.copy(current_interact_type))
                new_interacts_by_number[current_interact_type] = (atom_name, atom.mass_class)
                current_interact_type += 1

                match = new_interacts_by_mass_class[atom.mass_class]
                i = copy.copy(match[-1])
                j = copy.copy(match[-1])
                epsilon = opls_parameters.nonbonded[index].epsilon
                interact_type = InteractType(i, j, epsilon)
                new_interacts.append(interact_type)

            # next, add the charge types
            description = "\"OPLS %s (%s)\"" % (residue.name, atom_name)
            charge_type = ChargeType(current_charge_type, description, atom.charge)
            new_charges.append(charge_type)

            # now the biotypes
            biotype = BioType(current_biotype, atom_name, description, new_atom_types[idx].number, current_charge_type, current_bonded_type)
            residues_biotypes.append(biotype)
            new_biotypes.append(biotype)
            current_charge_type += 1
            current_biotype += 1
            current_bonded_type += 1

    # now, write the updated lines
    objects_to_write    = [new_atom_types, new_contacts, new_interacts, 
                           new_charges, new_bonds, new_angles, 
                           new_dihedrals, new_biotypes]

    absinth_locations   = [absinth_atom_locations, absinth_contact_locations, absinth_interact_locations, 
                           absinth_charge_locations, absinth_bond_locations, absinth_angle_locations,
                           absinth_torsion_locations, absinth_biotype_locations]
    names               = ['atom', 'contact', 'interact',
                           'charge', 'bonds', 'angles',
                           'torsion', 'biotype']
    header = '########## Added by OPLS-Parser ##########'
    start = 0
    savefile = open(savename, 'w')
    for objects, locations, name in zip(objects_to_write, absinth_locations, names):
        chunk = absinth_raw[start:locations[-1]+1] + [header] + output_objects_to_string(objects, name)
        start = locations[-1]+1
        savefile.write('\n'.join(chunk))
    savefile.close()