import os
import numpy
from tempfile import NamedTemporaryFile
from subprocess import check_output
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

def inclination(a, b, aerr, berr):
    i = numpy.arccos((1 + a**2 + b**2) ** (-1/2))
    ierr = numpy.sqrt(
        (aerr**2 * a + berr**2 * b) /
        ((1 - 1/(1 + a**2 + b**2)) * (1 + a**2 + b**2)**3)
    )
    return i, ierr

def position_angle(a, b, aerr, berr):
    theta = numpy.arctan(-a/b) + numpy.sign(b)*numpy.pi/2
    thetaerr = numpy.sqrt((aerr**2 / b**2 + berr**2 * a**2) /
                          ((a/b)**2 + 1)**2)

    return theta, thetaerr
