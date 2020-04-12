# Note: since OPLS actually uses the CHARMM format, it's not well-documented internally.
#
# This page was instrumental in understanding what the various terms mean: 
# https://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node25.html
#
# As was this one: http://cbio.bmt.tue.nl/pumma/index.php/Theory/Potentials
from common_objects import Mass, Atom
from common_objects import Bond, Angle, Dihedral, Improper, NonBonded, HBond
from common_objects import AtomType, BioType, BondLength, ChargeType, FOS, BondAngle, DihedralAngle
from common_objects import ContactType, InteractType, Radius
from common_objects import BondedTypeBond, BondedTypeTorsion, BondedTypeAngle, BondedTypeImptors
from common_objects import OPLS_Params, Absinth_Params
from common_objects import Residue, InternalCoordinate
from utils import strip_comments
from collections import OrderedDict
from pprint import pprint
from io import StringIO


# ---------------------------------------------------------------------------------------------------------------------

class Base_Parser(object):
    def __init__(self, filename):
        self._raw = list()
        self._data = list()
        self._comments_positions = list()

        self._read_file(filename)

        
    def _read_file(self, filename):
        data, data_no_comments, comments_positions = strip_comments(filename)
        self._raw = data
        self._data = data_no_comments
        self._comments_positions = comments_positions


    def _get_chunks(self, indices, collect_contiguous=False):
        # If this condition is true, we want to remove all the indices that are adjacent
        # as the current algorithm will create reduced copies.
        visit_indices = indices[::]
        if collect_contiguous:
            visit_indices = list()
            visit_indices.append(indices[0])
            for index, next_index in zip(indices[:-1], indices[1:]):
                diff_index = next_index - index

                # This is in case there is another match that is in another section of the file
                if diff_index != 1:
                    visit_indices.append(index)
        chunks = list()
        prev_index = -1
        for index in visit_indices:
            chunk = list()
            for line in self._data[index:]:
                if len(line.strip()) == 0:
                    break
                chunk.append(line)
            chunks.append(chunk)
        return chunks


    def _get_contiguous_chunks(self, indices):
        chunks = list()
        for index in indices:
            chunks.append(self._data[index])
        return chunks


    def parse(self):
        raise NotImplementedError('This method should not be called directly. It is meant to be overriden in a child class.')


# ---------------------------------------------------------------------------------------------------------------------

class OPLS_Topology_Parser(Base_Parser):
    def __init__(self, filename):
        super().__init__(filename)


    def _assign_masses(self, chunks):
        masses = OrderedDict()
        # There's only one mass section, hence the use of the 0 index.
        for chunk in chunks[0]:
            obj = chunk.split()
            number = int(obj[1])
            name = obj[2]
            mass = float(obj[3])
            element = obj[4]

            # ensure that only one record is added
            if name not in masses:
                masses[name] = Mass(number, name, mass, element)
            else:
                raise RuntimeError('Duplicate entries found for mass object "{name}"'.format(name=name))
        return masses


    def _parse_groups(self, groups_chunk):
        groups = list()
        temp = list()
        # group chunks always start with GROUP at the first index
        for line in groups_chunk[1:]:
            if line.startswith('ATOM'):
                ln = line.split()
                name, mass_class, charge = ln[1:]  # expand the slice into variables
                atom = Atom(name, mass_class, float(charge))
                temp.append(atom)
            else:
                groups.append(temp)
                temp = list()
        else:
            groups.append(temp)  # to get the last group of atoms (which is not terminated by GROUP)
        return groups


    def _parse_bond_angle_dihedral_chunk(self, chunk, connection, atoms_per_connection, error_message):
        all_connections = list()
        for line in chunk:
            if line.startswith(connection):
                ln = line.split()
                connections_in_line = ln[1:]  # the first token is the connection - e.g. "BOND", "ANGLE"
                assert len(connections_in_line) % atoms_per_connection == 0, 'The number of expected connections on a line must be divisible by {}.'.format(atoms_per_connection)

                connections = [connections_in_line[i*atoms_per_connection: (i*atoms_per_connection)+atoms_per_connection] 
                               for i in range(len(connections_in_line)//atoms_per_connection)]
                all_connections += connections
            else:
                raise RuntimeError(error_message)
        return all_connections


    def _parse_bonds(self, bonds_chunk):
        # Bonds are always pairs of atoms that are connected to each other
        all_bonds = list()
        for line in bonds_chunk:
            if line.startswith('BOND') or line.startswith('DOUBLE'):
                ln = line.split()
                bonds_in_line = ln[1:]  # the first token is "BOND" / "DOUBLE"
                assert len(bonds_in_line) % 2 == 0, 'The number of bonds on a line must be divisible by 2.'

                # Bonds occur in pairs. We can extract them by iterating over slices
                bonds = [bonds_in_line[i*2: (i*2)+2] for i in range(len(bonds_in_line)//2)]
                all_bonds += bonds
            else:
                raise RuntimeError('Invalid Bonds chunk - each line must start with BOND / DOUBLE.')
        return all_bonds


    def _parse_angles(self, angles_chunk):
        # Angles are 3 connected atoms
        error_message = 'Invalid ANGLE chunk - each line must start with ANGLE.'
        num_atoms_per_bond = 3
        all_angles = self._parse_bond_angle_dihedral_chunk(angles_chunk, 'ANGLE', num_atoms_per_bond, error_message)
        return all_angles


    def _parse_dihedrals(self, dihedrals_chunk):
        # Dihedrals are 4 atoms that are linearly connected
        error_message = 'Invalid DIHE chunk - each line must start with DIHE.'
        num_atoms_per_bond = 4
        all_dihedrals = self._parse_bond_angle_dihedral_chunk(dihedrals_chunk, 'DIHE', num_atoms_per_bond, error_message)
        return all_dihedrals


    def _parse_impr_bonds(self, impr_chunk):
        # Improper bond interactions are sets of 4 atoms for which an improper bond interaction will be calculated.
        # NOTE: The chiral centers is listed first.
        error_message = 'Invalid IMPR chunk - each line must start with IMPR.'
        num_atoms_per_bond = 4
        all_impr_bonds = self._parse_bond_angle_dihedral_chunk(impr_chunk, 'IMPR', num_atoms_per_bond, error_message)
        return all_impr_bonds


    def _parse_internal_coordinates(self, internal_coordinates_chunk):
        # Internal coordinates are those relative to other atoms in the residue. If atoms are missing from a structure
        # e.g. missing H-atoms in a PDB file - the complete molecule may still be built based on the positions of the
        # existing atoms.
        coordinates = list()
        for line in internal_coordinates_chunk:
            if not line.startswith('IC'):
                raise RuntimeError('Invalid IC chunk - each line must start with "IC".')
            ic = InternalCoordinate(line)
            coordinates.append(ic)
        return coordinates


    # -----------------------------------------------------------------------------------------------------------------
    # EXTRACTION SECTION

    def _extract_residue_information(self, chunk):
        name = None
        charge = 0.0
        is_patched = False
        for line in chunk:
            if line.startswith('RES') or line.startswith('PRES'):
                ln = line.split()[1:]
                name    = ln[0]
                charge  = float(ln[1])
                if line.startswith('PRES'):
                    is_patched = True
                break  # short-circuit (we don't need to parse the remainder of the chunk)
        return name, charge, is_patched


    def _extract_improper_bond_interactions(self, chunk):
        impr_chunk, impr_indices = self._extract_contiguous_parameter(chunk, 'IMPR')
        return impr_chunk, impr_indices


    def _extract_angles(self, chunk):
        angles_chunk, angles_indices = self._extract_contiguous_parameter(chunk, 'ANGLE')
        return angles_chunk, angles_indices


    def _extract_dihedrals(self, chunk):
        dihedrals_chunk, dihedrals_indices = self._extract_contiguous_parameter(chunk, 'DIHE')
        return dihedrals_chunk, dihedrals_indices


    def _extract_groups(self, chunk):
        groups_end_index = -1
        for index, line in enumerate(chunk):
            # find the first location of BOND or DELETE as GROUP always ends when these begin
            if line.startswith('BOND') or line.startswith('DELETE'):
                groups_end_index = index
                break

        # the first index of chunk is always information about the residue, patched or not
        groups_chunk = chunk[1:groups_end_index]
        groups_indices = list(range(1, groups_end_index))  # groups are enumerated in succession / contiguously
        return groups_chunk, groups_indices


    def _extract_contiguous_parameter(self, chunk, parameter):
        indices = list()
        trunc_chunk = list()
        for index, line in enumerate(chunk):
            if line.startswith(parameter):
                trunc_chunk.append(line)
                indices.append(index)
        return trunc_chunk, indices


    # Only PRES (patched residues have the DELETE line on them).
    def _extract_deleted_atoms(self, chunk):
        deleted_atoms_chunk, deleted_atoms_indices = self._extract_contiguous_parameter(chunk, 'DELETE')
        return deleted_atoms_chunk, deleted_atoms_indices


    def _extract_bonds(self, chunk):
        bonds_chunk, bonds_indices = self._extract_contiguous_parameter(chunk, 'BOND')
        double_bonds_chunk, double_bonds_indices = self._extract_contiguous_parameter(chunk, 'DOUBLE')  # double is a synonym for BOND
        combined_bonds_chunk = bonds_chunk + double_bonds_chunk
        combined_bonds_indices = bonds_indices + double_bonds_indices
        return combined_bonds_chunk, combined_bonds_indices


    def _extract_internal_coordinates(self, chunk):
        coordinates_chunk, coordinates_indices = self._extract_contiguous_parameter(chunk, 'IC')
        return coordinates_chunk, coordinates_indices


    def _extract_patches(self, chunk):
        patches_chunk, patches_indices = self._extract_contiguous_parameter(chunk, 'PATCHING')
        return patches_chunk, patches_indices


    def _extract_residues(self, chunks):
        '''Note: 
            BONDS are 2 connected atoms
            ANGLES are 3 connected atoms
            DIHEDRALS are 4 connected atoms (linearly)'''
        residues = OrderedDict()

        # first, break up the sections within each chunk into managable portions
        for chunk in chunks:
            name, charge, is_patched                    = self._extract_residue_information(chunk)
            groups_chunk, groups_indices                = self._extract_groups(chunk)
            deleted_atoms_chunk, deleted_atoms_indices  = self._extract_deleted_atoms(chunk)

            bonds_chunk, bonds_indices              = self._extract_bonds(chunk)
            angles_chunk, angles_indices            = self._extract_angles(chunk)
            dihedrals_chunk, dihedrals_indices      = self._extract_dihedrals(chunk)
            impr_chunk, impr_indices                = self._extract_improper_bond_interactions(chunk)
            coordinates_chunk, coordinates_indices  = self._extract_internal_coordinates(chunk)
            patches_chunk, patches_indices          = self._extract_patches(chunk)

            groups                  = self._parse_groups(groups_chunk)
            bonds                   = self._parse_bonds(bonds_chunk)
            impr_bonds              = self._parse_impr_bonds(impr_chunk)
            angles                  = self._parse_angles(angles_chunk)
            dihedrals               = self._parse_dihedrals(dihedrals_chunk)
            internal_coordinates    = self._parse_internal_coordinates(coordinates_chunk)

            residue = Residue(name, charge)
            residue.is_patched_residue = is_patched

            # Assign the atom groups and atoms to the residue
            for group_number, group in enumerate(groups, start=1):
                atom_groups = list()
                for atom in group:
                    atom_groups.append(atom.name)
                    residue.atoms[atom.name] = atom
                residue.atom_groups[group_number] = atom_groups
            
            # Populate the residue fields
            residue.bonds           = bonds[::]
            residue.improper_bonds  = impr_bonds[::]
            residue.angles          = angles[::]
            residue.dihedrals       = dihedrals[::]
            residue.coordinates     = internal_coordinates[::]

            # If the residue has a patch preference, parse it and assign
            if len(patches_chunk) == 1:
                residue.is_patched_residue  = True
                patches                     = patches_chunk[0].split()
                first_patch                 = patches[2]
                last_patch                  = patches[-1]
                if first_patch != 'NONE':
                    residue.patch_first = first_patch
                if last_patch != 'NONE':
                    residue.patch_last = last_patch
            
            # Add the residue
            residues[name] = residue
        return residues


    # This is where all the action happens
    def parse(self):
        decl_indices    = list()
        mass_indices    = list()
        res_indices     = list()
        pres_indices    = list()

        for index, line in enumerate(self._data):
            if line.startswith('DECL'): decl_indices.append(index)
            if line.startswith('MASS'): mass_indices.append(index)
            if line.startswith('RES'):  res_indices.append(index)
            if line.startswith('PRES'): pres_indices.append(index)

        decl_chunks = self._get_chunks(decl_indices)  # decl is one group / contiguous
        mass_chunks = self._get_chunks(mass_indices)  # mass is also one group / contiguous
        res_chunks = self._get_chunks(res_indices)
        pres_chunks = self._get_chunks(pres_indices)

        masses = self._assign_masses(mass_chunks)
        residues = self._extract_residues(res_chunks)
        patched_residues = self._extract_residues(pres_chunks)
        return masses, residues, patched_residues


# ---------------------------------------------------------------------------------------------------------------------


class OPLS_Parameters_Parser(Base_Parser):
    def __init__(self, filename):
        super().__init__(filename)

    def _get_chunks(self):
        # bond, angle, dihedral, improper, nonbonded, hbond
        temp = list()
        chunks = list()
        sections = 'BOND,ANGLE,DIHEDRAL,IMPROPER,NONBONDED,HBOND'.split(',')
        for line in self._data:
            ln = line.split()
            if len(ln) > 0:
                section = ln[0]
                if section not in sections:
                    temp.append(ln)
                else:
                    if len(temp) > 0:
                        chunks.append(temp[::])
                        temp = list()
                    if section == 'HBOND':
                        temp.append(ln)
                        chunks.append(temp)
                        temp = list()
        return chunks


    def _parse_bonds_chunk(self, bonds_chunk):
        all_bonds = list()
        for line in bonds_chunk:
            # recall that bonds are connections of two atoms!
            atom_a                  = line[0]
            atom_b                  = line[1]
            spring_constant         = float(line[2])
            reference_bond_length   = float(line[3])
            
            bond = Bond(atom_a, atom_b, spring_constant, reference_bond_length)
            all_bonds.append(bond)
        return all_bonds


    def _parse_angles_chunk(self, angles_chunk):
        all_angles = list()
        for line in angles_chunk:
            # recall that angles are connections of three atoms!
            atom_a                  = line[0]
            atom_b                  = line[1]
            atom_c                  = line[2]
            theta_force_constant    = float(line[3])
            reference_angle         = float(line[4])

            angle = Angle(atom_a, atom_b, atom_c, theta_force_constant, reference_angle)
            all_angles.append(angle)
        return all_angles


    def _parse_dihedrals_chunk(self, dihedrals_chunk):
        all_dihedrals = list()
        for line in dihedrals_chunk:
            # recall that dihedrals are connections of four atoms!
            atom_a              = line[0]
            atom_b              = line[1]
            atom_c              = line[2]
            atom_d              = line[3]
            chi_force_constant  = float(line[4])
            multiplicity        = int(line[5])
            delta               = float(line[6])

            dihedral = Dihedral(atom_a, atom_b, atom_c, atom_d, chi_force_constant, multiplicity, delta)
            all_dihedrals.append(dihedral)
        return all_dihedrals


    def _parse_improper_chunk(self, improper_chunk):
        all_improper = list()
        for line in improper_chunk:
            # recall that improper are connections of four atoms!
            atom_a              = line[0]
            atom_b              = line[1]
            atom_c              = line[2]
            atom_d              = line[3]
            chi_force_constant  = float(line[4])
            multiplicity        = int(line[5])
            delta               = float(line[6])

            improper = Improper(atom_a, atom_b, atom_c, atom_d, chi_force_constant, multiplicity, delta)
            all_improper.append(improper)
        return all_improper


    def _parse_nonbonded_chunk(self, nonbonded_chunk):
        all_nonbonded = list()
        # skip first line
        for line in nonbonded_chunk[1:]:
            atom        = line[0]
            ignored1    = float(line[1])
            epsilon     = float(line[2])
            rmin2       = float(line[3])
            ignored2    = float(line[4])
            epsilon_14  = float(line[5])
            rmin2_14    = float(line[6])

            nonbonded = NonBonded(atom, ignored1, epsilon, rmin2, ignored2, epsilon_14, rmin2_14)
            all_nonbonded.append(nonbonded)
        return all_nonbonded


    def _parse_hbonds_chunk(self, hbonds_chunk):
        cutoff = float(hbonds_chunk[0][-1])
        return HBond(cutoff)


    def parse(self):
        chunks = self._get_chunks()

        bonds_chunk     = chunks[0][::]
        angles_chunk    = chunks[1][::]
        dihedrals_chunk = chunks[2][::]
        improper_chunk  = chunks[3][::]
        nonbonded_chunk = chunks[4][::]
        hbonds_chunk    = chunks[5][::]

        bonds       = self._parse_bonds_chunk(bonds_chunk)
        angles      = self._parse_angles_chunk(angles_chunk)
        dihedrals   = self._parse_dihedrals_chunk(dihedrals_chunk)
        impropers   = self._parse_improper_chunk(improper_chunk)
        nonbonded   = self._parse_nonbonded_chunk(nonbonded_chunk)
        hbond       = self._parse_hbonds_chunk(hbonds_chunk)

        return OPLS_Params(bonds, angles, dihedrals, impropers, nonbonded, hbond)


# ---------------------------------------------------------------------------------------------------------------------


class Absinth_Parameters_Parser(Base_Parser):
    def __init__(self, filename):
        super().__init__(filename)
        self.angles                 = list()
        self.atoms                  = list()
        self.biotypes               = list()
        self.bonds                  = list()
        self.bonded_type_angles     = list()
        self.bonded_type_bonds      = list()
        self.bonded_type_imptors    = list()
        self.bonded_type_torsions   = list()
        self.charges                = list()
        self.contacts               = list()
        self.fos                    = list()
        self.interacts              = list()
        self.radii                  = list()
        self.torsions               = list()


    def _parse_angles(self, angles_chunk):
        for line in angles_chunk:
            ln = line.split()[1:]
            identifier = int(ln[0])
            potential_type = int(ln[1])
            a = float(ln[2])
            b = float(ln[3])
            c = 0.0
            d = 0.0
            if len(ln) == 6:
                c = float(ln[4])
                d = float(ln[5])
            angle = BondAngle(identifier, potential_type, a, b, c, d)
            self.angles.append(angle)


    def _parse_atoms(self, atoms_chunk):
        for line in atoms_chunk:
            ln = line.split()[1:]
            number = int(ln[0])
            name = ln[1]

            atomic_number = int(ln[-3])
            atomic_mass = float(ln[-2])
            valence = int(ln[-1])

            # remove the last elements ahead of joining for the description
            for i in range(3):
                ln.pop()

            description = ' '.join(ln[2:])[1:-1]
            atom = AtomType(number, name, description, atomic_number, atomic_mass, valence)
            self.atoms.append(atom)


    def _parse_biotypes(self, biotypes_chunk):
        # number, indicator, description, LJ_type, charge_type, bonded_type
        for line in biotypes_chunk:
            ln          = line.split()[1:]
            number      = int(ln[0])
            indicator   = ln[1]
            LJ_type     = int(ln[-3])
            charge_type = int(ln[-2])
            bonded_type = int(ln[-1])

            # remove the last elements ahead of joining for the description
            for i in range(3):
                ln.pop()

            description = ' '.join(ln[2:])[1:-1]
            biotype = BioType(number, indicator, description, LJ_type, charge_type, bonded_type)
            self.biotypes.append(biotype)


    def _parse_bonds(self, bonds_chunk):
        # identifier, potential_type, a, b, c
        for line in bonds_chunk:
            ln              = line.split()[1:]
            identifier      = int(ln[0])
            potential_type  = int(ln[1])
            a               = float(ln[2])
            b               = float(ln[3])
            c               = 0.0
            if len(ln) == 5:
                c = float(ln[4])  # not implemented?
            
            bond = BondLength(identifier, potential_type, a, b, c)
            self.bonds.append(bond)


    def _parse_bonded_type_angles(self, bonded_type_angles_chunk):
        for line in bonded_type_angles_chunk:
            ln = line.split()[1:]
            j               = int(ln[0])
            k               = int(ln[1])
            l               = int(ln[2])
            bond_potential  = int(ln[3])

            bond_angle = BondedTypeAngle(j, k, l, bond_potential)
            self.bonded_type_angles.append(bond_angle)


    def _parse_bonded_type_bonds(self, bonded_type_bonds_chunk):
        for line in bonded_type_bonds_chunk:
            ln = line.split()[1:]
            j               = int(ln[0])
            k               = int(ln[1])
            bond_potential  = int(ln[2])
            bond = BondedTypeBond(j, k, bond_potential)
            self.bonded_type_bonds.append(bond)


    def _parse_bonded_type_imptors(self, bonded_type_imptors_chunk):
        for line in bonded_type_imptors_chunk:
            ln = line.split()[1:]
            j               = int(ln[0])
            k               = int(ln[1])
            l               = int(ln[2])
            m               = int(ln[3])
            bond_potential  = int(ln[4])
            imptor = BondedTypeImptors(j, k, l, m, bond_potential)
            self.bonded_type_imptors.append(imptor)


    def _parse_bonded_type_torsions(self, bonded_type_torsions_chunk):
        for line in bonded_type_torsions_chunk:
            ln = line.split()[1:]
            j               = int(ln[0])
            k               = int(ln[1])
            l               = int(ln[2])
            m               = int(ln[3])
            bond_potential  = int(ln[4])
            torsion = BondedTypeTorsion(j, k, l, m, bond_potential)
            self.bonded_type_torsions.append(torsion)


    def _parse_charges(self, charges_chunk):
        # number, description, charge
        for line in charges_chunk:
            ln = line.split()[1:]
            number = int(ln[0])
            charge = float(ln[-1])
            ln.pop()
            description = ' '.join(ln[1:])[1:-1]
            charge = ChargeType(number, description, charge)
            self.charges.append(charge)


    def _parse_contacts(self, contacts_chunk):
        for line in contacts_chunk:
            ln = line.split()[1:]
            i           = int(ln[0])
            j           = int(ln[1])
            sigma_ij    = float(ln[2])
            contact = ContactType(i, j, sigma_ij)
            self.contacts.append(contact)


    def _parse_FOS(self, fos_chunk):
        for line in fos_chunk:
            ln = line.split()[1:]
            identifier      = ln[0]
            energy          = float(ln[1])
            enthalpy        = 0.0
            heat_capacity   = 0.0
            if len(ln) > 2:
                enthalpy = float(ln[2])
                heat_capacity = float(ln[3])
            fos = FOS(identifier, energy, enthalpy, heat_capacity)
            self.fos.append(fos)


    def _parse_interacts(self, interacts_chunk):
        for line in interacts_chunk:
            ln = line.split()[1:]
            i           = int(ln[0])
            j           = int(ln[1])
            epsilon_ij  = float(ln[2])
            interact = InteractType(i, j, epsilon_ij)
            self.interacts.append(interact)


    def _parse_radii(self, radii_chunk):
        for line in radii_chunk:
            ln = line.split()[1:]
            LJ_atom_type        = int(ln[0])
            replacement_radius  = float(ln[1])
            radius = Radius(LJ_atom_type, replacement_radius)
            self.radii.append(radius)


    def _parse_torsions(self, torsions_chunk):
        # identifier, potential_type, a, b, c, d, e, f, g
        for line in torsions_chunk:
            ln = line.split()[1:]
            identifier      = int(ln[0])
            potential_type  = int(ln[1])
            a = float(ln[2])
            b = float(ln[3])
            c = 0.0
            d = 0.0
            e = 0.0
            f = 0.0
            g = 0.0
            if len(ln) >= 8:
                c = float(ln[4])
                d = float(ln[5])
                e = float(ln[6])
                f = float(ln[7])
                if len(ln) == 9:
                    g = float(ln[8])
            torsion = DihedralAngle(identifier, potential_type, a, b, c, d, e, f, g)
            self.torsions.append(torsion)


    def parse(self):
        angle_indices               = list()
        atom_indices                = list()
        biotype_indices             = list()
        bond_indices                = list()
        bonded_type_angle_indices   = list()
        bonded_type_bond_indices    = list()
        bonded_type_imptors_indices = list()
        bonded_type_torsion_indices = list()
        charge_indices              = list()
        contact_indices             = list()
        fos_indices                 = list()
        interact_indices            = list()
        radius_indices              = list()
        torsion_indices             = list()
        for index, line in enumerate(self._data):
            ln = line.split()
            if len(ln) > 0:
                if ln[0] == 'angle':                angle_indices.append(index)
                if ln[0] == 'atom':                 atom_indices.append(index)
                if ln[0] == 'biotype':              biotype_indices.append(index)
                if ln[0] == 'bond':                 bond_indices.append(index)
                if ln[0] == 'bonded_type_angle':    bonded_type_angle_indices.append(index)
                if ln[0] == 'bonded_type_bond':     bonded_type_bond_indices.append(index)
                if ln[0] == 'bonded_type_imptors':  bonded_type_imptors_indices.append(index)
                if ln[0] == 'bonded_type_torsion':  bonded_type_torsion_indices.append(index)
                if ln[0] == 'charge':               charge_indices.append(index)
                if ln[0] == 'contact':              contact_indices.append(index)
                if ln[0] == 'fos':                  fos_indices.append(index)
                if ln[0] == 'interact':             interact_indices.append(index)
                if ln[0] == 'radius':               radius_indices.append(index)
                if ln[0] == 'torsion':              torsion_indices.append(index)

        angles_chunk                = self._get_contiguous_chunks(angle_indices)
        atoms_chunk                 = self._get_contiguous_chunks(atom_indices)
        biotypes_chunk              = self._get_contiguous_chunks(biotype_indices)
        bonds_chunk                 = self._get_contiguous_chunks(bond_indices)
        bonded_type_angles_chunk    = self._get_contiguous_chunks(bonded_type_angle_indices)
        bonded_type_bonds_chunk     = self._get_contiguous_chunks(bonded_type_bond_indices)
        bonded_type_imptors_chunk   = self._get_contiguous_chunks(bonded_type_imptors_indices)
        bonded_type_torsions_chunk  = self._get_contiguous_chunks(bonded_type_torsion_indices)
        charges_chunk               = self._get_contiguous_chunks(charge_indices)
        contacts_chunk              = self._get_contiguous_chunks(contact_indices)
        fos_chunk                   = self._get_contiguous_chunks(fos_indices)
        interacts_chunk             = self._get_contiguous_chunks(interact_indices)
        radii_chunk                 = self._get_contiguous_chunks(radius_indices)
        torsions_chunk              = self._get_contiguous_chunks(torsion_indices)

        self._parse_angles(angles_chunk)
        self._parse_atoms(atoms_chunk)
        self._parse_biotypes(biotypes_chunk)
        self._parse_bonds(bonds_chunk)
        self._parse_bonded_type_angles(bonded_type_angles_chunk)
        self._parse_bonded_type_bonds(bonded_type_bonds_chunk)
        self._parse_bonded_type_imptors(bonded_type_imptors_chunk)
        self._parse_bonded_type_torsions(bonded_type_torsions_chunk)
        self._parse_charges(charges_chunk)
        self._parse_contacts(contacts_chunk)
        self._parse_FOS(fos_chunk)
        self._parse_interacts(interacts_chunk)
        self._parse_radii(radii_chunk)
        self._parse_torsions(torsions_chunk)

