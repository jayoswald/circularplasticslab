#!/usr/bin/env python3
import numpy
from numpy import pi


class lammps_system:

    def __init__(self):
        self.atom_coord = []
        # Molecule ID for each atom.
        self.mol_id = []
        self.atom_type = []
        # List of tuples defining each bond: (bond type, i1, i2).
        self.bonds = []
        self.box = []
        self.masses = []


    def num_atoms(self):
        ''' Returns the number of atoms in the system '''
        return len(self.atom_type)

    def build_bond_table(self):
        ''' Builds table of atoms bonded to each atom for faster bond lookups. '''
        bt = [[] for _ in range(num_atoms())]
        for (i1, i2, _) in self.bonds:
            bt[i1].append(i2)
            bt[i2].append(i1)
        for row in bt:
            row.sort()
        return bt

    def write_lmpdata(self, path):
        ''' Writes the system to a LAMMPS data file. '''
        fid = open(path, 'w')
        fid.write('LAMMPS polymer data file\n\n')
        fid.write(f'{self.num_atoms():>12d}  atoms\n')
        fid.write(f'{len(self.bonds):>12d}  bonds\n\n')
        atom_types = set(self.atom_type)
        bond_types = {b[0] for b in self.bonds}
        fid.write(f'{len(atom_types):>12d}  atom types\n')
        fid.write(f'{len(bond_types):>12d}  bond types\n\n')
        fid.write(f'{self.box[0]:>15f} {self.box[1]:>15f} xlo xhi\n')
        fid.write(f'{self.box[2]:>15f} {self.box[3]:>15f} ylo yhi\n')
        fid.write(f'{self.box[4]:>15f} {self.box[5]:>15f} zlo zhi\n')
        fid.write('\nMasses\n\n')
        for i, m in enumerate(self.masses):
            fid.write(f'{i+1:>12d} {m:9f}\n')
        fid.write('\nAtoms\n\n')
        for i, x in enumerate(self.atom_coord):
            fid.write(f'{i+1:>7d} {self.mol_id[i]:6d} {self.atom_type[i]:6d} ')
            fid.write(f'{x[0]:9.6f} {x[1]:9.6f} {x[2]:9.6f}\n')
        fid.write('\nBonds\n\n')
        for i, (t, j, k) in enumerate(self.bonds):
            fid.write(f'{i+1:>7d} {t:>6d} {j+1:>9d} {k+1:>9d}\n')



def construct_random_chains(**kwargs):
    ''' '''
    # Number of molecules in the system.
    molecules = kwargs.get('molecules', 10)
    # How many beads per molecule.
    molecule_size = kwargs.get('molecule_size', 20)

    n = molecules*molecule_size
    v0 = 4.0/3.0*pi
    Lx = (n / v0)**(1.0/3.0)
    Lx = kwargs.get('box_length_x', Lx)
    Ly = kwargs.get('box_length_y', Lx)
    Lz = kwargs.get('box_length_z', Lx)
    
    bond_length = kwargs.get('bond_length', 1.0)
    system = lammps_system()

    system.box = [0.0, Lx, 0.0, Ly, 0.0, Lz]
    system.masses = [1]

    for mol in range(molecules):
        for i in range(molecule_size):
            if i == 0:
                # First atom in chain
                x = numpy.random.rand(3) * [Lx, Ly, Lz]
            elif i == 1:
                # Second atom in chain
                x = bond_length*random_unit() + system.atom_coord[-1]
            else:
                t = system.atom_coord[-1] - system.atom_coord[-2]
                t /= numpy.linalg.norm(t)
                n = random_perpendicular(t)
                # TODO - sample angle q
                R = rotation_about_axis(n, 2.0*pi/3.0)
                x = system.atom_coord[-1] + R@t*bond_length
            
            system.atom_coord.append(x)
            system.atom_type.append(1)
            system.mol_id.append(mol)
            if i > 0:
                system.bonds.append((1, system.num_atoms()-2, system.num_atoms()-1))

    return system

def random_unit():
    ''' Returns a unit vector with random orientation. '''
    x = numpy.random.normal(size=3)
    return x / numpy.linalg.norm(x)


def random_perpendicular(v):
    ''' Returns a random unit vector perpendicular to v. '''
    v /= numpy.linalg.norm(v)
    while True:
        u = random_unit()
        if abs(numpy.dot(u, v)) < 0.9:
            break
    p = numpy.cross(u, v)
    return p / numpy.linalg.norm(p)


def rotation_about_axis(u, q):
    ''' Returns a rotation matrix for a rotation of q (radians) about axis u. '''
    x, y, z = u / numpy.linalg.norm(u)
    c = numpy.cos(q)
    s = numpy.sin(q)
    R = [[c+x*x*(1.0-c),   x*y*(1.0-c)-z*s, x*z*(1.0-c)+y*s],
         [y*x*(1.0-c)+z*s, c+y*y*(1.0-c),   y*z*(1.0-c)-x*s],
         [z*x*(1.0-c)-y*s, z*y*(1.0-c)+x*s, c+z*z*(1.0-c)]]
    return numpy.array(R)


if __name__ == '__main__':
    system = construct_random_chains()
    system.write_lmpdata('test.lammps')
