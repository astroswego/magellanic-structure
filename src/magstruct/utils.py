import numpy

def demean(X):
    """Subtract the mean from each column of X"""
    for i in range(X.shape[1]):
        X[:,i] = X[:,i] - X[:,i].mean()
    return X

def arctan2(a, b):
    return numpy.arctan(a/b) + (0          if a > 0 and b > 0 else
                                numpy.pi   if b < 0 else
                                2*numpy.pi)
