#!/usr/bin/env python

import numpy as np
import sys

if __name__ == '__main__':
    """Calculates sum of .cube voxels with value >= threshold. 

       Usage:  cube_sum_greater.py density.cube threshold"""

    threshold = float(sys.argv[2])

    with open(sys.argv[1]) as f:
        f.readline()
        f.readline()
        nidxp = int(f.readline().split()[0])
        nx = int(f.readline().split()[0])
        ny = int(f.readline().split()[0])
        nz = int(f.readline().split()[0])
        for i in range(nidxp):
            f.readline()
        a = np.fromfile(f, sep=' ').reshape(nx,ny,nz)
        #b = a[np.nonzero(a)]

        sum_gt = a[a>=threshold].sum()
        print "{} of {} = {}".format(sum_gt, a.sum(), sum_gt/a.sum())
