from argparse import ArgumentError, ArgumentParser, FileType
from magstruct import plane, transformations
from numpy import array, float as npfloat, loadtxt, savetxt
from os import path
import matplotlib.pyplot as plt

def get_args():
    parser = ArgumentParser()
    parser.add_argument('-i', '--input', type=FileType('r'),
        help='file to read RA/Dec/Distance values from')
    parser.add_argument('-o', '--output', type=str,
        help='directory to output plots to')
    parser.add_argument('-c', '--center', type=npfloat, nargs=3,
        metavar='d RA Dec',
        help='coordinates to center of target galaxy')

    args = parser.parse_args()

    return args

def main():
    args = get_args()

    ra_dec_dist = loadtxt(args.input)
    xyz = array([transformations.radec2xyz(rdd, args.center)
                 for rdd in ra_dec_dist])
    savetxt(path.join(args.output, "cartesian.dat"), xyz)
    x, y, z = xyz.T
    plt.scatter(x, z)
    plt.savefig(path.join(args.output, "cartesian.png"))
    plt.clf()
    ra, dec, dist = ra_dec_dist.T
    plt.scatter(ra, dec)
    plt.savefig(path.join(args.output, "radec.png"))

if __name__ == '__main__':
    exit(main())
