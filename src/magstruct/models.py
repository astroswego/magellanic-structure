import numpy
from sklearn.linear_model import LinearRegression

class Plane():
    def __init__(self, regressor=LinearRegression(fit_intercept=True)):
        self.regressor = regressor

    def fit(self, X, y=None):
        assert X.shape[1] == 3, 'X must have 3 columns: x, y, and z'

        self.regressor.fit(X[:,:2], X[:,2])
        return self

    def transform(self, X):
        assert X.shape[1] == 3, 'X must have 3 columns: x, y, and z'

        self.regressor.predict(X[:,:2])

        self.coef = numpy.empty(3)
        self.coef[:2] = self.regressor.coef_
        self.coef[2] = self.regressor.intercept_

        return self.coef, self.rotations(self.coef)

    def score(X, y=None, sample_weight=None):
        return self.regressor.score(X[:,:2], X[:,2], sample_weight)

    @staticmethod
    def rotations(coef):
        """Finds the inclination and position angles for a plane defined by
        the given coefficients.
        """
        a, b, c = coef

        i = numpy.arccos((1 + a**2 + b**2)**(-1/2))
        theta = numpy.nan

        return i, theta


class Ellipsoid():
    def __init__(self):
        pass

    def fit(self, X, y=None):
        return self

    def transform(self, X, y=None):
        self.I = self.inertia_tensor(X)

        evals, evecs = numpy.linalg.cov(self.I)


    @staticmethod
    def inertia_tensor(X):
        x, y, z = X.T

        cov_xx = numpy.cov(x, bias=0)
        cov_yy = numpy.cov(y, bias=0)
        cov_zz = numpy.cov(z, bias=0)
        cov_xy = numpy.cov(x, y, bias=0)[0][1]
        cov_xz = numpy.cov(x, z, bias=0)[0][1]
        cov_yz = numpy.cov(y, z, bias=0)[0][1]

        I_xx = cov_yy + cov_zz
        I_yy = cov_xx + cov_zz
        I_zz = cov_xx + cov_yy
        I_xy = I_yx = cov_xy
        I_xz = I_zx = cov_xz
        I_yz = I_zy = cov_yz

        return numpy.array([
            [ I_xx, -I_xy, -I_xz],
            [-I_yx,  I_yy,  I_yz],
            [-I_zx, -I_zy,  I_zz]
        ])
