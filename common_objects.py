from collections import OrderedDict, namedtuple

# ---------------------------------------------------------------------------------------------------------------------
# OPLS Topology objects

# No need for a class here when a named tuple can do
Mass = namedtuple('Mass', 'number, name, mass, element')    # for topology
Atom = namedtuple('Atom', 'name, mass_class, charge')       # for topology

# ---------------------------------------------------------------------------------------------------------------------
# OPLS Parameter objects

Bond        = namedtuple('Bond', 'atom_a, atom_b, spring_constant, reference_bond_length')
Angle       = namedtuple('Angle', 'atom_a, atom_b, atom_c, theta_force_constant, reference_angle')
Dihedral    = namedtuple('Dihedral', 'atom_a, atom_b, atom_c, atom_d, chi_force_constant, multiplicity, delta')
Improper    = namedtuple('Improper', 'atom_a, atom_b, atom_c, atom_d, psi_force_constant, multiplicity, delta')
NonBonded   = namedtuple('NonBonded', 'atom, ignored1, epsilon, Rmin2, ignored2, epsilon_14, Rmin2_14')
HBond       = namedtuple('HBond', 'cutoff')

# ---------------------------------------------------------------------------------------------------------------------
# ABSINTH objects

# See: http://campari.sourceforge.net/V3/parameters.html#S2_LJ-types
AtomType            = namedtuple('AtomType', 'number, name, description, atomic_number, atomic_mass, valence')

# See: http://campari.sourceforge.net/V3/parameters.html#S1_Biotypes
BioType             = namedtuple('BioType', 'number, indicator, description, LJ_type, charge_type, bonded_type')

# See: http://campari.sourceforge.net/V3/parameters.html#S5_Bond_length_potential_types
BondLength          = namedtuple('BondLength', 'identifier, potential_type, a, b, c')

# See: http://campari.sourceforge.net/V3/parameters.html#S3_Charge_types
ChargeType          = namedtuple('ChargeType', 'number, description, charge')

# See: http://campari.sourceforge.net/V3/parameters.html#S2_LJ-types
ContactType         = namedtuple('ContactType', 'i, j, sigma_ij')

# See: http://campari.sourceforge.net/V3/parameters.html#S4_Free_energies_of_solvation
FOS                 = namedtuple('FOS', 'identifier, energy, enthalpy, heat_capacity')

# See: http://campari.sourceforge.net/V3/parameters.html#S6_Bond_angle_potential_types
BondAngle           = namedtuple('BondAngle', 'identifier, potential_type, a, b, c, d')

# See: http://campari.sourceforge.net/V3/parameters.html#S7_Dihedral_angle_torsional_potential
DihedralAngle       = namedtuple('DihedralAngle', 'identifier, potential_type, a, b, c, d, e, f, g')

# See: http://campari.sourceforge.net/V3/parameters.html#S2_LJ-types
InteractType        = namedtuple('InteractType', 'i, j, epsilon_ij')

# See: http://campari.sourceforge.net/V3/parameters.html#S2_LJ-types
Radius              = namedtuple('Radius', 'LJ_atom_type, replacement_radius')

# See: http://campari.sourceforge.net/V3/parameters.html#S9_Bonded_types
BondedTypeBond      = namedtuple('BondedTypeBond', 'j, k, bond_potential')
BondedTypeAngle     = namedtuple('BondedTypeAngle', 'j, k, l, bond_potential')  # only jkl and lkj are recognized permutations
BondedTypeTorsion   = namedtuple('BondedTypeTorsion', 'j, k, l, m, bond_potential')  # potential from torsion
BondedTypeImptors   = namedtuple('BondedTypeImptors', 'j, k, l, m, bond_potential')  # potetial from torsion

# ---------------------------------------------------------------------------------------------------------------------
# Return objects

OPLS_Params = namedtuple('OPLS_Params', 'bonds, angles, dihedrals, impropers, nonbonded, hbond')
Absinth_Params = namedtuple('Absinth_Params', 'angles, atoms, '
                                              'biotypes, bonds, bonded_type_angles, bonded_type_bonds, bonded_type_imptors, '
                                              'bonded_type_torsions, charges, contacts, fos, interacts, radii, torsions')

# ---------------------------------------------------------------------------------------------------------------------

# This is a Plain-Ole-Data (POD) object
class InternalCoordinate(object):
    def __init__(self, coordinates_data):
        self.atom_A = None
        self.atom_B = None
        self.atom_C = None
        self.atom_D = None
        self.is_improper = False  # set if `atom_C` is improper (has an asterisk next to its name)

        self.bond_length_AB = 0.0
        self.bond_length_CD = 0.0
        self.bond_length_AC = 0.0  # only set if  `atom_C` is improper (alternates with AB)

        self.angle_ABC = 0.0
        self.angle_BCD = 0.0
        self.angle_BCA = 0.0  # only set if `atom_C` is improper (alternates with ABC)

        self.dihedral_ABCD = 0.0
        self._parse_coordinate_data(coordinates_data)


    def _parse_coordinate_data(self, coordinates_data):
        ln = coordinates_data.split()[1:]  # The first token is always "IC"
        self.atom_A = ln[0]
        self.atom_B = ln[1]
        self.atom_C = ln[2]
        self.atom_D = ln[3]
        if '*' in self.atom_C:
            self.is_improper = True
            self.bond_length_AC = float(ln[4])
            self.angle_BCA      = float(ln[5])
        else:   
            self.bond_length_AB = float(ln[4])
            self.angle_ABC      = float(ln[5])

        self.dihedral_ABCD  = float(ln[6])
        self.angle_BCD      = float(ln[7])
        self.bond_length_CD = float(ln[8])


    def __repr__(self):
        if self.is_improper:
            class_repr = 'InternalCoordinate(A={}, B={}, C={}, D={}, AC-length={}, BCA-angle={}, ABCD-dihedral={}, BCD-angle={}, CD-length={})'
            return class_repr.format(self.atom_A, self.atom_B, self.atom_C, self.atom_D,
                                     self.bond_length_AC, self.angle_BCA, 
                                     self.dihedral_ABCD, self.angle_BCD, self.bond_length_CD)
        else:
            class_repr = 'InternalCoordinate(A={}, B={}, C={}, D={}, AB-length={}, ABC-angle={}, ABCD-dihedral={}, BCD-angle={}, CD-length={})'
            return class_repr.format(self.atom_A, self.atom_B, self.atom_C, self.atom_D,
                                     self.bond_length_AB, self.angle_ABC, 
                                     self.dihedral_ABCD, self.angle_BCD, self.bond_length_CD)

    def __str__(self):
        return 'InternalCoordinate(A={}, B={}, C={}, D={}'.format(self.atom_A, self.atom_B, self.atom_C, self.atom_D)


# ---------------------------------------------------------------------------------------------------------------------


# A generic class meant to store attributes for both unpatched and patched residues.
class Residue(object):
    def __init__(self, name=None, charge=0.0):
        self.name = name

        # groups are collection of atom objects
        self.atom_groups = OrderedDict()
        self.bonds = list()
        self.improper_bonds = list()
        self.angles = list()
        self.dihedrals = list()

        # atoms are defined from mass objects but ascribed a common name and charge
        self.atoms = OrderedDict()
        self.deleted_atoms = OrderedDict()

        # the internal coordinates of the atoms
        self.coordinates = OrderedDict()

        # For patched residues
        self.is_patched_residue = False

        # This is to cover the case where the line "PATCHING FIRS <PRES> LAST <PRES>" appears at the end (for both RES and PRES)
        self.has_patching_preference = False
        self.patch_first = None
        self.patch_last = None
        self.charge = charge


    def _find_atom(self, atom_name):
        name = atom_name[::]
        modifier = atom_name[0]
        if modifier in '+-*':
            name = atom_name[1:]
        return self.atoms[name]


    def _assign_atoms_from_container(self, connections_container):
        # A connection can be a bond, angle, or dihedral
        atomized = OrderedDict()
        for connection in connections_container:
            user_friendly_name = ' - '.join(connection)
            atomized[user_friendly_name] = tuple(self._find_atom(atom_name) for atom_name in connection)
        return atomized


    # This and the other `get` methods were intended to be used to populate atoms. However, as-is, it cannot handle missing atoms
    # TODO: Fix!
    def get_bonds(self):
        '''This returns an OrderedDict in which the atom names within the bonds are populated by Atom objects.'''
        return self._assign_atoms_from_container(self.bonds)


    def get_angles(self):
        '''This returns an OrderedDict in which the atom names within the angles are populated by Atom objects.'''
        return self._assign_atoms_from_container(self.angles)


    def get_dihedrals(self):
        '''This returns an OrderedDict in which the atom names within the dihedrals are populated by Atom objects.'''
        return self._assign_atoms_from_container(self.dihedrals)


    def __repr__(self):
        residue = 'Residue(name={name}, charge={charge}, num-atoms={num_atoms}, atom-groups={num_groups}, num-bonds={num_bonds}, angles={angles}, dihedrals={dihedrals}, patched={patched})'
        return residue.format(name=self.name, 
                              charge=self.charge,
                              num_atoms=len(self.atoms),
                              num_groups=len(self.atom_groups),
                              num_bonds=len(self.bonds),
                              angles=len(self.angles),
                              dihedrals=len(self.dihedrals),
                              patched=self.is_patched_residue)


    def __str__(self):
        return 'Residue({name}, {charge})'.format(name=self.name, charge=self.charge)
