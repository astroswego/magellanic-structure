from argparse import ArgumentError, ArgumentParser, FileType
from magstruct import plane, plots, transformations
from numpy import array, float as npfloat, loadtxt, savetxt
from os import path
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import pi, meshgrid, dot
from magstruct.utils import demean

rad2deg = 180/pi

def get_args():
    parser = ArgumentParser()
    parser.add_argument('-i', '--input', type=FileType('r'),
        help='file to read RA/Dec/Distance values from')
    parser.add_argument('-o', '--output', type=str,
        help='directory to output plots to')
    parser.add_argument('-c', '--center', type=npfloat, nargs=3,
        metavar=('RA', 'DEC', 'DIST'),
        help='coordinates to center of target galaxy')
    parser.add_argument('-r', '--regression', type=str, default='LTS',
        choices=['OLS', 'LTS'],
        help='regression method to use')
    parser.add_argument('--xlim', type=float, default=None, metavar='X',
        help='limit on distance from center in x of stars in plane fitting')
    parser.add_argument('--ylim', type=float, default=None, metavar='Y',
        help='limit on distance from center in y of stars in plane fitting')
    parser.add_argument('--zlim', type=float, default=None, metavar='Z',
        help='limit on distance from center in z of stars in plane fitting')

    args = parser.parse_args()

    args.fit = plane.fit if args.regression == 'OLS' else plane.lts_fit

    return args

def main():
    args = get_args()

    ra_dec_dist = loadtxt(args.input)
    xyz = array([transformations.radec2xyz(rdd, args.center)
                 for rdd in ra_dec_dist])
    demean(xyz) # subtract the mean from each column
    savetxt(path.join(args.output, "cartesian.dat"), xyz)
    x, y, z = xyz.T
    plots.plot2d(x, y, path.join(args.output, "unrotated_xy.png"),
                 title='Unrotated cartesian')

    ra, dec, dist = ra_dec_dist.T
    plots.plot2d(ra, dec, path.join(args.output, "radec.png"),
                 xlabel='RA', ylabel='Dec', title='RA/Dec')

    (a, b, c), (aerr, berr, cerr) = args.fit(xyz, **args.__dict__)
    print("a = {} +- {}, b = {} +- {}, c = {} +- {}".format(
        a, aerr, b, berr, c, cerr))
    plots.plot3d(x, y, z, a, b, c)
    i, ierr = plane.inclination(a, b, aerr, berr)
    theta, thetaerr = plane.angle(a, b, aerr, berr)
    print("i = {} +- {}, theta = {} +- {}".format(i*rad2deg, ierr*rad2deg,
                                                  theta*rad2deg,
                                                  thetaerr*rad2deg))

    xyz_ = array([transformations.rotate2plane(x, y, z, i, theta)
                  for x, y, z in xyz])
    x_, y_, z_ = xyz_.T
    plots.plot2d(x_, y_, path.join(args.output, "rotated_xy.png"),
                 title='Rotated cartesian')
    plots.plot2d(x_, z_, path.join(args.output, "rotated_xz.png"),
                 title='Rotated cartesian')

if __name__ == '__main__':
    exit(main())
