import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

__all__ = [
    'plot3d'
]

def plot2d(x, y, filename,
           alpha=1.0, color=[0,0,1],
           xlabel='x', ylabel='y', title=''):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(x, y, alpha=alpha, color=color)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    fig.savefig(filename)
    fig.clf()
           
def plot3d(x, y, z, a, b, c,
           scatter_alpha=1.0, scatter_color=[0,0,1],
           surface_alpha=0.2, surface_color=[1,0,0],
           xlabel='x', ylabel='y', zlabel='z', title=''):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, z,
               alpha=scatter_alpha, color=scatter_color)
    point = np.array([0.0, 0.0, c])
    normal = np.array(np.cross([1, 0, a], [0, 1, b]))
    d = -point.dot(normal)
    xs, ys = np.meshgrid(np.arange(x.min(), x.max(), 0.1),
                         np.arange(y.min(), y.max(), 0.1))
    zs = (-normal[0] * xs - normal[1] * ys - d) / normal[2]
    ax.plot_surface(xs, ys, zs,
                    alpha=surface_alpha, color=surface_color)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    ax.set_title(title)
    plt.show()
    fig.clf()
