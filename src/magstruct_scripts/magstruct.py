from argparse import ArgumentError, ArgumentParser, FileType, SUPPRESS
from magstruct import models, transformations
from numpy import array, float as npfloat, loadtxt, savetxt
from os import path
import matplotlib.pyplot as plt
import numpy
from sklearn.pipeline import Pipeline

def get_args():
    parser = ArgumentParser()
    parser.add_argument('-i', '--input', type=FileType('r'),
        help='file to read RA/Dec/Distance values from')
    parser.add_argument('-o', '--output', type=str,
        help='directory to output plots to')
    parser.add_argument('-c', '--center', type=npfloat, nargs=3,
        metavar=('RA', 'DEC', 'DIST'),
        help='coordinates to center of target galaxy')
    parser.add_argument('-r', '--regression', type=str, default='OLS',
        choices=['OLS'],
        help='regression method to use')
    parser.add_argument('--xlim', type=float, default=None, metavar='X',
        help='limit on distance from center in x of stars in plane fitting')
    parser.add_argument('--ylim', type=float, default=None, metavar='Y',
        help='limit on distance from center in y of stars in plane fitting')
    parser.add_argument('--zlim', type=float, default=None, metavar='Z',
        help='limit on distance from center in z of stars in plane fitting')

    args = parser.parse_args()

    return args

def main():
    args = get_args()

    equatorial_coordinates = loadtxt(args.input)
    RA_0, Dec_0, D_0 = args.center

    # transforms coordinates from equatorial to cartesian
    eq2cart = transformations.Equatorial2Cartesian(RA_0, Dec_0, D_0)

    # plane fitting pipeline
    plane = Pipeline([('Cartesian', eq2cart),
                      ('Plane', models.Plane())])
    # ellipsoid fitting pipeline
    ellipsoid = Pipeline([('Cartesian', eq2cart),
                          ('Ellipsoid', models.Ellipsoid())])

    (a, b, c), (i, theta) = plane.fit_transform(equatorial_coordinates)

    print('(a, b, c) =', ', '.join(map(str,(a, b, c))))
    print('(i, theta) =', ', '.join(map(lambda x: str(numpy.rad2deg(x)),
                                        (i, theta))))

    
    exit()

    
    cartesian_coordinates = eq2cart.transform(equatorial_coordinates)

#    demean(xyz) # subtract the mean from each column
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
#    plots.plot3d(x, y, z, a, b, c)
    i, ierr = plane.inclination(a, b, aerr, berr)
    theta, thetaerr = plane.position_angle(a, b, aerr, berr)
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
    I = ellipsoid.moment_of_inertia_tensor(xyz)
    unsorted_eigvals, unsorted_eigvecs = eig(I)
    index_array = numpy.fromiter(
        reversed(unsorted_eigvals.argsort(kind='mergesort')),
        dtype=int
    )
    eigvals = unsorted_eigvals[index_array]
    eigvecs = unsorted_eigvecs[index_array,:]
    exit(print(eigvals))
    S_1, S_2, S_3 = ellipsoid.get_axes(eigvals)
    i, ierr = ellipsoid.inclination(eigvecs)
    theta, thetaerr = ellipsoid.position_angle(eigvecs)
    print("S_1 = {}, S_2 = {}, S_3 = {}".format(S_1, S_2, S_3))
    print("i = {} +- {}, theta = {} +- {}".format(i*rad2deg, None,#ierr*rad2deg,
                                                  theta*rad2deg,
                                                  None))#thetaerr*rad2deg))

if __name__ == '__main__':
    exit(main())
