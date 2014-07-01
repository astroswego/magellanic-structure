import numpy
from sklearn.linear_model import LinearRegression

def fit(xyz):
    xy = xyz[:,:2]
    z  = xyz[:,-1]
    ones = numpy.ones((len(z), 1))
    xy1 = numpy.hstack((xy, ones))
    ols = LinearRegression()
    ols.fit(xy1, z)
    coeffs = ols.coef_
    coeffs[-1] = ols.intercept_
    fit = ols.predict(xy1)

    return coeffs, fit

def inclination(a, b):
    return numpy.arccos((1 + a**2 + b**2) ** (-1/2))

def angle(a, b):
    return numpy.arctan(-a/b) + numpy.sign(b)*numpy.pi/2
