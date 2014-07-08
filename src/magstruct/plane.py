import numpy
from statsmodels.api import add_constant, OLS
from lts_fits.planefit import PlaneFit

def fit(xyz, xlim=None, ylim=None, zlim=None, **kwargs):    
    all_true = numpy.empty_like(xyz[:,0], dtype=bool) \
               if None in [xlim, ylim, zlim] \
               else None
    xbool = numpy.abs(xyz[:,0]) < xlim if xlim else all_true
    ybool = numpy.abs(xyz[:,1]) < ylim if ylim else all_true
    zbool = numpy.abs(xyz[:,2]) < zlim if zlim else all_true
    bools = numpy.logical_and(numpy.logical_and(xbool, ybool), zbool)
    XYZ = xyz[bools,:]
    XY = add_constant(XYZ[:,:2], prepend=False)
    Z  = XYZ[:,-1]
    model = OLS(Z, XY)
    result = model.fit()
    coeffs = result.params
    stderr = result.HC1_se

    return coeffs, stderr

def lts_fit(xyz, xlim=None, ylim=None, zlim=None, **kwargs):
    all_true = numpy.empty_like(xyz[:,0], dtype=bool) \
               if None in [xlim, ylim, zlim] \
               else None
    xbool = numpy.abs(xyz[:,0]) < xlim if xlim else all_true
    ybool = numpy.abs(xyz[:,1]) < ylim if ylim else all_true
    zbool = numpy.abs(xyz[:,2]) < zlim if zlim else all_true
    bools = numpy.logical_and(numpy.logical_and(xbool, ybool), zbool)
    x, y, z = xyz[bools,:].T
    # eventually put real errors here
    xnoise = ynoise = znoise = numpy.random.normal(0, 0.1, x.size)
    lts_fit = PlaneFit(x, y, z,
                       xnoise, ynoise, znoise,
                       pivotx=numpy.median(x), pivoty=numpy.median(y))
    c, a, b = lts.abc
    cerr, aerr, berr = lts.abc_err

    return [a,b,c], [aerr, berr, cerr]

def inclination(a, b, aerr, berr):
    i = numpy.arccos((1 + a**2 + b**2) ** (-1/2))
    ierr = numpy.sqrt(
        (aerr**2 * a + berr**2 * b) /
        ((1 - 1/(1 + a**2 + b**2)) * (1 + a**2 + b**2)**3)
    )
    return i, ierr

def angle(a, b, aerr, berr):
    theta = numpy.arctan(-a/b) + numpy.sign(b)*numpy.pi/2
    thetaerr = numpy.sqrt((aerr**2 / b**2 + berr**2 * a**2) /
                          ((a/b)**2 + 1)**2)

    return theta, thetaerr
