#!/usr/bin/env python

import argparse
import warnings
import MDAnalysis as mda
import MDAnalysis.topology.core as topcore
from MDAnalysis.core.topologyattrs import Bonds, Angles, Dihedrals, Impropers


def find_impropers_by_type(u, t):
    """Find impropers with atom types (t0, t1, t2, t3).
       Selects angles of type (t0, t2, t3) having other bond (t2,t1)."""
    dihedrals_found = set()
    for a in u.angles.select_bonds((t[0], t[2], t[3])):
        atom = a[1]  # select middle atom in angle
        for other_bond in atom.bonds:
            other_atom = other_bond.partner(atom)
            # if this atom is of right type and isn't in the angle we started with
            if other_atom.type == t[1] and other_atom not in a:
                # if searching with a[1], want tuple of (a[0], new, a[1], a[2])
                # check that angle was not reversed
                if a[0].type == t[0]:
                    desc = (a[0].index, other_atom.index, a[1].index, a[2].index)
                else:
                    desc = (a[2].index, other_atom.index, a[1].index, a[0].index)
                if desc[0] > desc[-1]:
                    desc = desc[::-1]
                dihedrals_found.add(desc)

    return tuple(dihedrals_found)


def build_zif8_top(fname, guess_all_impropers=False):
    """Build periodic ZIF-8 topology from PDB file."""

    mda.core.flags['use_pbc'] = True
    u = mda.Universe(fname, guess_bonds=False)

    # set atom types from PDB atom names for proper bond typing
    u.atoms.types = u.atoms.names

    # guess bonds using periodic box
    rvdw = {'Zn': 1.4, 'N': 1.5, 'C1': 1.4, 'C2': 1.4, 'C3': 1.4, 'H2': 1.0, 'H3': 1.0}
    all_bonds = topcore.guess_bonds(u.atoms, u.atoms.positions, vdwradii=rvdw, box=u.dimensions)
    u.add_TopologyAttr(Bonds(all_bonds))

    # guess angles
    all_angles = topcore.guess_angles(u.atoms.bonds)
    u.add_TopologyAttr(Angles(all_angles))

    # guess dihedrals
    all_dihedrals = topcore.guess_dihedrals(u.atoms.angles)
    u.add_TopologyAttr(Dihedrals(all_dihedrals))

    # guess impropers
    if guess_all_impropers:
        impropers = topcore.guess_improper_dihedrals(u.atoms.angles)
    else:
        # MDanalysis does not guarantee the order of atoms in the result of guess_improper_dihedrals(), 
        # so by now we use custom search of predefined dihedral types
        # select Zheng's impropers
        ia1 = find_impropers_by_type(u, ('N', 'C3', 'C1', 'N'))
        ia2 = find_impropers_by_type(u, ('C2', 'H2', 'C2', 'N'))
        ia3 = find_impropers_by_type(u, ('C2', 'Zn', 'N', 'C1'))
        impropers = ia1 + ia2 + ia3
    u.add_TopologyAttr(Impropers(impropers))

    return u


def canon_bondstr(t):
    """Return canonical string representaion for bond, angle, etc.
    It has the form of AtomName1-AtomName2-...
    Use reverse order if it is lexicographically smaller."""
    if tuple(t[::-1]) < tuple(t):
        t = t[::-1]

    return '-'.join(t)


def print_atoms(u):
    print '[ moleculetype ]'
    print '; Name       nrexcl'
    print 'ZIF    3'
    print
    print '; %d atoms of %d types' % (len(u.atoms), len(set(u.atoms.types)))
    print '[ atoms ]'
    print ';   nr       type  resnr residue  atom   cgnr     charge       mass'
    for a in u.atoms:
        print "%6d %10s %5d %6s %6s %7d" % (a.id, a.name + '_zif8', a.resnum, 'ZIF', a.name, a.id)
    print


def print_bonds(u):
    print '; %d bonds of %d types' % (len(u.bonds), len(u.bonds.types()))
    print ';', ', '.join([canon_bondstr(t) for t in u.bonds.types()])
    print '[ bonds ]'
    print ';  ai    aj  funct    b0        kb'
    for b in u.bonds:
        btype = 'gb_zif8_' + canon_bondstr(b.type)
        print "%5s %5s      %-20s ; %-4s %-4s" % (b[0].id, b[1].id, btype, b[0].name, b[1].name)
    print


def print_pairs(u):
    # find 1-4 pairs
    pairs_all = set([tuple(sorted([i1, i4])) for i1, i2, i3, i4 in u.dihedrals.to_indices()])
    exclusions = set([tuple(sorted([i1, i3])) for i1, i2, i3 in u.angles.to_indices()])
    pairs14 = pairs_all - exclusions

    print '; %d pairs' % (len(pairs14))
    print '[ pairs ]'
    print ';  ai    aj   funct'
    for i1, i2 in sorted(pairs14):
        print "%5d %5d %5d    ; %-3s %-3s" % (i1+1, i2+1, 1, u.atoms[i1].name, u.atoms[i2].name)
    print


def print_angles(u):
    print '; %d angles of %d types' % (len(u.angles), len(u.angles.types()))
    print ';', ', '.join([canon_bondstr(t) for t in u.angles.types()])
    print '[ angles ]'
    print ';  ai    aj    ak   funct      phi0      fc'
    for a in u.angles:
        atype = 'ga_zif8_' + canon_bondstr(a.type)
        print "%5s %5s %5s       %-20s ; %-4s %-4s %-4s" % (a[0].id, a[1].id, a[2].id, atype, a[0].name, a[1].name, a[2].name)
    print


def print_dihedrals(u):
    print '; %d diherals of %d types' % (len(u.dihedrals), len(u.dihedrals.types()))
    print ';', ', '.join([canon_bondstr(t) for t in u.dihedrals.types()])
    print '[ dihedrals ]'
    print ';  ai    aj    ak    al   funct   phi0        fc    mult'
    for d in u.dihedrals:
        dtype = 'gd_zif8_' + canon_bondstr(d.type)
        print "%5s %5s %5s %5s     %-20s ; %-3s %-3s %-3s %-3s" % (d[0].id, d[1].id, d[2].id, d[3].id,
                                                                   dtype, d[0].name, d[1].name, d[2].name, d[3].name)
    print


def print_impropers(u):
    print '; %d impropers of %d types' % (len(u.impropers), len(u.impropers.types()))
    print ';', ', '.join([canon_bondstr(t) for t in u.impropers.types()])
    print '[ dihedrals ]'
    print ';  ai    aj    ak    al   funct    phi0         fc'
    for d in u.impropers:
        dtype = 'gi_zif8_' + canon_bondstr(d.type)
        print "%5s %5s %5s %5s    %-20s ; %-3s %-3s %-3s %-3s" % (d[0].id, d[1].id, d[2].id, d[3].id,
                                                                  dtype, d[0].name, d[1].name, d[2].name, d[3].name)
    print


if __name__ == '__main__':

    tested_MDAnalysis_version = '0.17.0'
    parser = argparse.ArgumentParser(description='Generates GROMACS topology include file from ZIF-8 PDB structure file.')
    parser.add_argument('fname', default='conf.pdb', metavar='conf.pdb', nargs='?', help='ZIF-8 PDB file')
    parser.add_argument('-a', '--all-impropers', dest='all_impropers', action='store_true', help='guess all improper dihedrals')
    args = parser.parse_args()

    if mda.version.__version__ != tested_MDAnalysis_version:
        warnings.warn("This scripts was only tested with MDAnalysis version {}.".format(tested_MDAnalysis_version))

    u = build_zif8_top(args.fname, guess_all_impropers=args.all_impropers)

    print '; Generated by pdb2top_ZIF-8.py script from {} file.\n'.format(args.fname)
    print_atoms(u)
    print_bonds(u)
    print_pairs(u)
    print_angles(u)
    print_dihedrals(u)
    print_impropers(u)
