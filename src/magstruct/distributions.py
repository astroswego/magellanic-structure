from itertools import product, combinations
import numpy
from numpy import sqrt
from numpy.random import choice, uniform

# def uniform_cube(step=0.01):
#     return array(list(product(arange(-1,1,step), repeat=3)))

def uniform_ellipsoid(a, b, c, i=0, j=0, k=0, step=0.01):
    """Generates a grid of points which fit within an ellipsoid"""
    limit = max(a,b,c)
    return array([[x,y,z]
                  for x,y,z in product(arange(-limit, limit, step), repeat=3)
                  if ((x-i)/a)**2 + ((y-j)/b)**2 + ((z-k)/c)**2 <= 1
    ])

def random_ellipsoid(a, b, c, i=0, j=0, k=0, scatter=(0.5, 1), step=0.01):
    uni_ellipsoid = uniform_ellipsoid(a, b, c, i, j, k, step)
    noise = uniform(scatter[0], scatter[1], uni_ellipsoid.shape)
    pm = choice([-1,1], size=noise.shape, replace=True)

    return uni_ellipsoid + pm*noise

# def _in_ellipsoid_fn(a,b,c):
#     return lambda xyz: (
#         xyz[0]**2/a**2 + xyz[1]**2/b**2 + xyz[2]**2/c**2 <= 1
#     )

def ellipsoid_distribution(a, b, c, n=1000, wmax=None):
    assert a**2 > b**2 > c**2

    if wmax is None:
        wmax = -c**2 + (c**2)/2
    else:
        assert wmax > -c**2

    u = uniform(-b**2, -c**2, n)
    v = uniform(-a**2, -b**2, n)
    w = uniform(-c**2,  wmax, n)

    return _uvw2x(a,b,c,u,v,w), _uvw2y(a,b,c,u,v,w), _uvw2z(a,b,c,u,v,w)

def _uvw2x(a, b, c, u, v, w):
    pm = choice([-1,1], size=u.shape, replace=True)
    return pm*sqrt((a**2+w)*(a**2+u)*(a**2+v)
                   /(a**2-b**2)/(a**2-c**2))

def _uvw2y(a, b, c, u, v, w):
    pm = choice([-1,1], size=u.shape, replace=True)
    return pm*sqrt((b**2+w)*(b**2+u)*(b**2+v)
                   /(b**2-a**2)/(b**2-c**2))

def _uvw2z(a, b, c, u, v, w):
    pm = choice([-1,1], size=u.shape, replace=True)
    return pm*sqrt((c**2+w)*(c**2+u)*(c**2+v)
                   /(c**2-b**2)/(c**2-a**2))

def is_ellipsoidal(a, b, c):
    return lambda X: (
        (X[:,0]/a)**2 + (X[:,1]/b)**2 + (X[:,2]/c)**2 <= 1
    )

def monte_carlo(f, shape, low=0.0, high=1.0, scale=3/2, dtype=float, order='C'):
    assert scale > 1.0
    samples = numpy.empty(shape, dtype=dtype)
    nsamples = 0
    required, variables = shape
    while nsamples < required:
        new_generated = numpy.multiply(numpy.subtract(shape, [nsamples, 0]),
                                       [scale, 1])
        generated = numpy.random.uniform(low, high, new_generated)
        satisfied = generated[f(generated), :]
        needed = required-nsamples
        new = min(needed, satisfied.shape[0])
        samples[nsamples:nsamples+new, :] = satisfied[:new, :]
        nsamples += new
    return samples

def monte_carlo_ellipsoid(a, b, c, N, *args, **kwargs):
    f = is_ellipsoidal(a,b,c)
    shape = (N, 3)
    semi_major_axis = max(a,b,c)
    return monte_carlo(f, shape, low=-semi_major_axis, high=semi_major_axis,
                       *args, **kwargs)
