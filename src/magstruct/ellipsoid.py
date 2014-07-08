import numpy

def moment_of_inertia_tensor(xyz):
    x, y, z = xyz.T

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

def diagonalize(I):
    eigvals, eigvecs = numpy.linalg.eig(I)

    return eigvals * numpy.eye(3)
