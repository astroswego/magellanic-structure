import numpy

def demean(X):
    """Subtract the mean from each column of X"""
    for i in range(X.shape[1]):
        X[:,i] = X[:,i] - X[:,i].mean()
    return X
