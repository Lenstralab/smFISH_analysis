import numpy as np
from sklearn import mixture
from tllab_common.misc import *

def gaussian_unmix(d, n):
    """ get the components in a gaussian mixture
        d: data, 1D
        n: number of components

        A, mu, s: arrays with amplitude, mu and sigma for each component
    """
    m = mixture.GaussianMixture(n)
    m.fit(d.reshape(-1, 1))
    A = m.weights_ / np.sqrt(2*np.pi * m.covariances_.flatten())
    mu = m.means_.flatten()
    s = np.sqrt(m.covariances_).flatten()
    return A, mu, s


def gaussian(x, a, m, s):
    return a * np.exp(-(((x-m)/s)**2)/2)


def outliers(D, keep=True):
    """ get the indices of outliers in D
        keep=True: get the indices of everything but the outliers
    """
    D = np.array(D).flatten()
    q2 = np.nanmedian(D)
    q1 = np.nanmedian(D[D < q2])
    q3 = np.nanmedian(D[D > q2])
    lb = 4 * q1 - 3 * q3
    ub = 4 * q3 - 3 * q1
    if keep:
        idx = np.where((D >= lb) & (D <= ub))[0]
    else:
        idx = np.where(~((D >= lb) & (D <= ub)))[0]
    return idx