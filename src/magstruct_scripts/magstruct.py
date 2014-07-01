from argparse import ArgumentError, ArgumentParser, FileType
from magstruct import plane, transformations
from numpy import array, float as npfloat, loadtxt, savetxt
from os import path
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def get_args():
    parser = ArgumentParser()
    parser.add_argument('-i', '--input', type=FileType('r'),
        help='file to read RA/Dec/Distance values from')
    parser.add_argument('-o', '--output', type=str,
        help='directory to output plots to')
    parser.add_argument('-c', '--center', type=npfloat, nargs=3,
        metavar=('RA', 'DEC', 'DIST'),
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
    plt.scatter(x, y)
    plt.savefig(path.join(args.output, "unrotated.png"))
    plt.clf()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, z)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()

    ra, dec, dist = ra_dec_dist.T
    plt.scatter(ra, dec)
    plt.savefig(path.join(args.output, "radec.png"))
    plt.clf()

    (a, b, c), fitted_values = plane.fit(xyz)
    i, theta = plane.inclination(a, b), plane.angle(a, b)

    xyz_ = array([transformations.rotate2plane(x, y, z, i, theta)
                  for x, y, z in xyz])
    x_, y_, z_ = xyz_.T
    plt.scatter(x_, y_)
    plt.savefig(path.join(args.output, "rotated_xy.png"))
    plt.clf()

    plt.scatter(x_, z_)
    plt.savefig(path.join(args.output, "rotated_xz.png"))
    plt.clf()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x_, y_, z_)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()

if __name__ == '__main__':
    exit(main())
