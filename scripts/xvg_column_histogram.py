#!/usr/bin/env python

import numpy as np
import argparse


if __name__ == '__main__':
    """Compute histogram of selected column"""

    parser = argparse.ArgumentParser(description='Compute histogram of selected column of .xvg file', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', dest='fname_xvg', default='dist.xvg', metavar='dist.xvg', help='xvg file')
    parser.add_argument('-c', dest='col', type=int, default=2, metavar=2, help='column number')
    parser.add_argument('-o', dest='fname_out', default='hist.dat', metavar='hist.dat', help='output file with histogram')
    parser.add_argument('-nonorm', action='store_true', help='Do not normalize')
    parser.add_argument('-bins', type=int, default=50, help='Number of bins')
    parser.add_argument('-xmin', type=float, default=0, help='min range')
    parser.add_argument('-xmax', type=float, default=1.0, help='max range (if rmax=rmin use actual min/max values')

    args = parser.parse_args()

    x = np.loadtxt(args.fname_xvg, comments=('#','@'))
    print "Loaded .xvg data with shape: ", x.shape
    if len(x.shape) == 1:
        x = x.reshape(-1, 1)
        print " reshaped to: ", x.shape
    col = args.col - 1
    hist_data = x[:,col]
    data_min = hist_data.min()
    data_max = hist_data.max()

    print "hist data shape: ", hist_data.shape
    print "first, last values: {}, {}".format(hist_data[0], hist_data[-1])
    print "min, max data values: {}, {}".format(data_min, data_max)

    if args.xmin == args.xmax:
        args.xmin = data_min
        args.xmax = data_max

    hist, bins = np.histogram(hist_data, bins=args.bins, range=(args.xmin, args.xmax), density=not args.nonorm)
    bincenters = 0.5*(bins[1:]+bins[:-1])
    hdr = "mean: {}  std: {}".format(hist_data.mean(), hist_data.std())
    np.savetxt(args.fname_out, np.array((bincenters,hist)).T, fmt="%g", header=hdr, delimiter=' ')
