import numpy
from numpy import array, sin, cos

__all__ = [
    'Equatorial2Cartesian',
    'Rotation3D',
    'rotation_matrix_3d'
]

class Equatorial2Cartesian():
    def __init__(self, RA_0, Dec_0, D_0):
        self.RA_0 = RA_0
        self.Dec_0 = Dec_0
        self.D_0 = D_0

    def fit(self, X, y=None):
        return self

    def transform(self, X, y=None, **params):
        X_new = numpy.empty_like(X)

        x, y, z = X_new[:,0], X_new[:,1], X_new[:,2]
        RA, Dec, D = X[:,0], X[:,1], X[:,2]

        delta_RA = RA - self.RA_0

        x[:] = -D * sin(delta_RA) * cos(Dec)
        y[:] =  D * (sin(Dec) * cos(self.Dec_0) +
                     sin(self.Dec_0) * cos(delta_RA) * cos(Dec))
        z[:] = self.D_0 \
             - D * (sin(Dec)*sin(self.Dec_0) + cos(RA)*cos(self.Dec_0)) \
             - self.RA_0*cos(Dec)

        return X_new

def rotation_matrix_3d(angle, axis):
    assert axis in range(3), 'Axis must be 0, 1, or 2'

    T = numpy.empty((3, 3), dtype=float)

    # find the index of the -sin(angle) term
    # this formula is the polynomial which passes through all of the pairs
    # (axis, index)
    i = axis**2 - 4*axis + 5
    
    T.flat[::3+1] = cos(angle)
    T.flat[i::3-1] = sin(angle)
    # negate the -sin(angle) term, as it is currently just sin(angle)
    T.flat[i] *= -1
    
    T[axis,:] = 0
    T[:,axis] = 0
    T[axis,axis] = 1

    return T
    
class Rotation3D():
    def __init__(self, angle, axis):
        self.axis = axis
        self.angle = angle
        self.rotation_matrix = rotation_matrix_3d(angle, axis)

    def fit(self, X, y=None):
        return self

    def transform(self, X, y=None, **params):
        return self.rotation_matrix.dot(X)
