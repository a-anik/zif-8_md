#!/usr/bin/env python

import argparse
import numpy as np

# ZIF-8 with all H, space group 197'
space_group = 'I 2 3'
a = 16.991
b = 16.991
c = 16.991

residues = [ ('LIG', [
('C1', np.array([ 0.517921, 0.124949 , 0.124949   ]) ),
('C2', np.array([ 0.404685, 0.132188 , 0.188335   ]) ),
('C2', np.array([ 0.404685, 0.188335 , 0.132188   ]) ),
('C3', np.array([ 0.595139, 0.0949326, 0.0949326  ]) ),
('N' , np.array([ 0.475252, 0.0921664, 0.185687   ]) ),
('N' , np.array([ 0.475252, 0.185687 , 0.0921664  ]) ),
('H2', np.array([ 0.359249, 0.120946 , 0.227003   ]) ),
('H2', np.array([ 0.359249, 0.227003 , 0.120946   ]) ),
('H3', np.array([ 0.600965, 0.109411 , 0.0339003  ]) ),
('H3', np.array([ 0.600965, 0.0339003, 0.109411   ]) ),
('H3', np.array([ 0.642046, 0.121535 , 0.12742    ]) ),
]),

('ZN', [
('Zn', np.array([ 0.5000, 0.2500,  0         ]) ),
])
]

# lattice vectors
e1 = np.array([a, 0, 0])
e2 = np.array([0, b, 0])
e3 = np.array([0, 0, c])
E = np.array([e1, e2, e3])


def gen_symcopies(res_atoms):

    # symmetry operations from .cif file
    op = dict()
    op[1] = lambda x,y,z: (x,y,z)
    op[2] = lambda x,y,z: (0.5+x,0.5+y,0.5+z)
    op[3] = lambda x,y,z: (-x,y,-z)
    op[4] = lambda x,y,z: (-x,-y,z)
    op[5] = lambda x,y,z: (x,-y,-z)
    op[6] = lambda x,y,z: (0.5-x,0.5+y,0.5-z)
    op[7] = lambda x,y,z: (0.5-x,0.5-y,0.5+z)
    op[8] = lambda x,y,z: (0.5+x,0.5-y,0.5-z)
    op[9] = lambda x,y,z: (z,x,y)
    op[10] = lambda x,y,z: (y,z,x)
    op[11] = lambda x,y,z: (-z,-x,y)
    op[12] = lambda x,y,z: (z,-x,-y)
    op[13] = lambda x,y,z: (-z,x,-y)
    op[14] = lambda x,y,z: (y,-z,-x)
    op[15] = lambda x,y,z: (-y,z,-x)
    op[16] = lambda x,y,z: (-y,-z,x)
    op[17] = lambda x,y,z: (0.5+z,0.5+x,0.5+y)
    op[18] = lambda x,y,z: (0.5-z,0.5-x,0.5+y)
    op[19] = lambda x,y,z: (0.5+z,0.5-x,0.5-y)
    op[20] = lambda x,y,z: (0.5-z,0.5+x,0.5-y)
    op[21] = lambda x,y,z: (0.5+y,0.5+z,0.5+x)
    op[22] = lambda x,y,z: (0.5+y,0.5-z,0.5-x)
    op[23] = lambda x,y,z: (0.5-y,0.5+z,0.5-x)
    op[24] = lambda x,y,z: (0.5-y,0.5-z,0.5+x)

    out = []
    for sym_op in [op[i] for i in range(1,25)]:
        out.append([(name, sym_op(*tuple(v))) for name, v in res_atoms])

    return out


def write_atom(serial=None, name=None, altLoc=None, resName=None, chainID=None,
               resSeq=None, iCode=None, x=None, y=None, z=None, occupancy=1.0, tempFactor=0.0,
               segID=None, element=None, charge=""):
    """Write PDB ATOM record.

    Adapted from MDAnalyis package (http://www.mdanalysis.org).
    """

    for arg in ('serial', 'name', 'resName', 'resSeq', 'x', 'y', 'z',
                'occupancy', 'tempFactor', 'charge'):
        if locals()[arg] is None:
            raise ValueError('parameter ' + arg + ' must be defined.')
    serial = int(str(serial)[-5:])  # check for overflow here?
    name = name[:4]
    if len(name) < 4:
        name = " " + name   # customary to start in column 14
    altLoc = altLoc or " "
    altLoc = altLoc[:1]
    resName = resName[:4]
    chainID = chainID or ""   # or should we provide a chainID such as 'A'?
    chainID = chainID.strip()[-1:]  # take the last character
    resSeq = int(str(resSeq)[-4:])  # check for overflow here?
    iCode = iCode or ""
    iCode = iCode[:1]
    element = element or aname.rstrip("0123456789")        # element == 0|False|None will be deduced from name
    element = str(element).strip()[:2]                     # make sure that is a string for user input
    segID = segID or chainID
    segID = segID[:4]
    fmt = "ATOM  %(serial)5d %(name)-4s %(resName)-4s%(chainID)1s%(resSeq)4d%(iCode)1s   %(x)8.3f%(y)8.3f%(z)8.3f%(occupancy)6.2f%(tempFactor)6.2f      %(segID)-4s%(element)2s%(charge)2s"
    print fmt % vars()


def remove_duplicates(res_copies):
    # coordinates of all residue copies: ncopies x natoms x 3
    x = np.array([np.array([acoord for aname, acoord in r]) for r in res_copies])
    # put coords to primary box
    x = np.remainder(x, 1.0)
    out = []
    for i, r in enumerate(res_copies):
        # compare residue copy with all previous and keep unique
        if not any([np.allclose(x[i], x[j], atol=0.001) for j in range(i)]):
            out.append(r)

    return out


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Build supercell of ZIF-8.')
    parser.add_argument('n1', type=int, help='number of cell X-repetitions')
    parser.add_argument('n2', type=int, help='number of cell Y-reperitions')
    parser.add_argument('n3', type=int, help='number of cell Z-repetitions')
    args = parser.parse_args()

    # number of elementary cells in each direction
    n1, n2, n3 = args.n1, args.n2, args.n3

    print "MODEL     1"
    print "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4i" % (a*n1, b*n2, c*n3, 90, 90, 90, space_group, 1)

    # init sequential atom number and residue number for writing
    atnr = 1
    resnr = 1
    supercell_shifts = [(e1*i + e2*j + e3*k) for i in range(n1) for j in range(n2) for k in range(n3)]

    for resname, resatoms in residues:
        sym_copies = gen_symcopies(resatoms)
        sym_copies = remove_duplicates(sym_copies)
        for t in supercell_shifts:
            for res in sym_copies:
                for aname, fract_coord in res:
                    v = np.dot(fract_coord, E)+t
                    write_atom(serial=atnr, name=aname, resName=resname, resSeq=resnr, x=v[0], y=v[1], z=v[2])
                    atnr += 1
                resnr += 1

    print "END"
