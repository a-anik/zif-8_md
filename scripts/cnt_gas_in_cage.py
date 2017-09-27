#!/usr/bin/env python

import numpy as np
import argparse
import MDAnalysis
from scipy.spatial import ConvexHull


def molecule_in_hull(r, A, b):
    """Test key atoms/points of the residue r if they are inside the convex polytope (A,b)"""

    #r.pack_into_box()   # pbc pack, no need with typical trajectories, slow
    if r.name == 'CO2':
        #positions = [r.O1.position, r.O2.position]
        positions = [r.C.position]
    elif r.name == 'N2':
        #positions = [r.N1.position, r.N2.position]
        positions = [r.M.position]

    #positions = [r.center_of_geometry()]   # r must be pbc-whole for meaningful cog calculation
    #positions = r.positions   # all atoms

    positions_in_hull = [np.all(np.dot(A, x) + b < 0) for x in positions]
    
    # Alternative is to return np.any(positions_in_hull)
    return np.all(positions_in_hull)


if __name__ == '__main__':
    """Counts CO2/N2 molecules in central cage of 2x2x2 ZIF-8 framework"""

    parser = argparse.ArgumentParser(description='Counts number of gas molecules in central cage of ZIF-8 MOF', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', dest='fntraj', default='traj.xtc', metavar='traj.xtc', help='trajectory')
    parser.add_argument('-s', dest='fntop', default='conf.gro', metavar='conf.gro', help='topology file')
    parser.add_argument('-o', dest='fnout', default='gas_mols_cage0.xvg', metavar='gas_mols_cage0.xvg', help='output file')
    parser.add_argument('-oh', dest='fnhistout', default='hist_gas_mols_cage0.dat', metavar='hist_gas_mols_cage0.dat', help='output histogram file')
    args = parser.parse_args()

    u = MDAnalysis.Universe(args.fntop, args.fntraj)

    # select Zn atoms forming central cage with TEMPO inside it
    bx, by, bz, _, _, _ = u.dimensions
    central_cage = u.select_atoms("name Zn and point {} {} {} {}".format(bx/2, by/2, bz/2, 11.0)) # Zn within 11A of the box center
    if len(central_cage) != 24:
        raise ValueError("Central cage does not have 24 vertices of truncated octahedron")

    # select gas molecules
    sel = u.select_atoms("resname CO2 or resname N2")

    hist_data = []
    fout = open(args.fnout, 'w')
    
    for ts in u.trajectory:
        ##import MDAnalysis.analysis.distances
        ##o = np.array([bx/2, by/2, bz/2], dtype=np.float32).reshape(1,3)
        ##box = np.array([bx,by,bz], dtype=np.float32)
        ##dists = MDAnalysis.analysis.distances.distance_array(o, central_cage.positions, box)

        # construct convex hull around cage vertices
        hull = ConvexHull(central_cage.positions)
        A = hull.equations[:,0:3]
        b = hull.equations[:,3]

        # count gas molecule inside the cage
        num_in_cage = np.sum([molecule_in_hull(mol, A, b) for mol in sel.residues])
        hist_data.append(num_in_cage)

        fout.write("{:.2f} {}\n".format(ts.time, num_in_cage))

    fout.close()

    hist_data = np.array(hist_data)
    max_val = hist_data.max()

    hist, bins = np.histogram(hist_data, bins=max_val+1, range=(-0.5, max_val+0.5), density=True)
    bincenters = 0.5*(bins[1:]+bins[:-1])
    hdr = "mean: {}  std: {}".format(hist_data.mean(), hist_data.std())
    np.savetxt(args.fnhistout, np.array((bincenters,hist)).T, fmt="%g", header=hdr, delimiter=' ')

