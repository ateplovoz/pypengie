'''
datatools.py

Module with various data and statistic-related functions

---
Created: Vadim Pribylov, 2018-10-19
'''

import numpy as np
from scipy.stats import chi2, t

BASEALPHA = 0.05
BASEZSTAR = 2.325


def zstar_help():
    '''
    Returns z* table
    '''
    return (
            '| alpha |   z*  |\n'
            '|-------|-------|\n'
            '|  99%  | 2.576 |\n'
            '|  95%  | 2.326 |\n'
            '|  95%  | 1.960 |\n'
            '|  90%  | 1.645 |\n'
        )


def clean3sigma(data, mask=True):
    '''
    Removes points with deviation from mean more than three sigma

    Args:
        data: numpy.ndarray - input data
        mask: bool - whether to return boolean mask

    Returns:
        numpy.ndarray - array or bool mask of cleaned data
    '''
    try:
        data_mean = np.mean(data)
        data_std = np.std(data)
    except TypeError:
        print('error: expected numpy.ndarray type, got {}'.format(type(data)))
        return np.nan
    data_mask = np.abs(data - data_mean) <= data_std
    if mask:
        return data_mask
    else:
        return data[data_mask]


def var_interv(data):
    '''
    Returns range midpoint
    '''
    try:
        midpoints = np.array([
            (data[i] + data[i + 1]) / 2 for i in range(len(data) - 1)
        ])
    except TypeError:
        print('error: expected list type, got {}'.format(type(data)))
        return np.nan
    return midpoints


def phi(self, u):
    '''
    phi(u) --- normal distribution parameter.

    Not sure what it does.

    Args:
        u - a parameter of some sort
    '''
    # TODO: fix the docstring
    return 1 / np.sqrt(2*np.pi) * np.power(np.e, -np.square(u)/2)


def chi2test(data, datatype='data', alpha=BASEALPHA):
    '''
    Pearson's chi-squared test (χ2) for normal distribution check

    Args:
        data: np.ndarray - input data
        datatype: str - type of input data. Can be `data` or `hist`.

    Returns:
        dict - dictionary of parameters, including histogram ranges
    '''
    if dtype == 'data':
        (freq, ranges) = np.histogram(data)
    elif dtype == 'hist':
        try:
            (freq, ranges) = data
        except TypeError:
            raise TypeError(
                'error: expected numpy.ndarray type, got {}'.format(type(data))
            )
        except ValueError:
            raise ValueError(
                'error: expected data to by tuple of (freq, ranges), '
                'got {0} of length of {1}'.format(type(data), len(data))
            )
    else:
        raise AttributeError('unsupported attribute')
    med_var = var_interv(ranges)
    h = med_var[1] - med_var[0]
    xav = np.average(med_var, weights=freq)
    # weighted STD
    xstd = np.sqrt(
        np.average(np.square(med_var), weights=freq)
        - np.square(xav)
    )
    ranges_i = np.sum(freq) * h / xstd * phi((med_var-xav) / xstd)
    chi2exp = np.sum(np.square(freq-ranges_i) / ranges_i)
    # order of freedem
    k = med_var.size - 3
    # chi-squared critical value
    chi2Cr = chi2.ppf(1 - alpha, k)
    return {
        'freq': freq,
        'ranges': ranges,
        'ranges_i': ranges_i,
        'chi2cr': chi2Cr,
        'chi2': chi2exp,
        'av': xav,
        'std': xstd,
        'H0': chi2exp < chi2Cr
        'H1': chi2exp > chi2Cr
    }


def conf_level(data, method='std'):
    '''
    Calculation of confidence level

    Args:
        data: np.ndarray - input data
        method: str - what method use for calculation. `var` for variance
        or `std` for standard deviation
    '''
    conf_lower = None
    conf_upper = None
    tstud = None
    try:
        mean = np.mean(data)
        std = np.std(data)
        variance = np.var(data)
    except TypeError:
        raise TypeError(
            'error: expected numpy.ndarray type, got {}'.format(type(data))
        )
    n = data.size
    if method == 'var':
        # calculate using variance
        tstud = t.ppf(1 - alpha/2, n - 1)
        conf_lower = mean - tstud*variance/np.sqrt(n)
        conf_upper = mean + tstud*variance/np.sqrt(n)
    elif method == 'std':
        # calculate using standard deviation
        tstud = t.ppf(1 - self.alpha/2, n - 1)
        conf_range = tstud * std / np.sqrt(n)
        conf_lower = mean - conf_range
        conf_upper = mean + conf_range
    else:
        raise AttributeError('unsupported `method`. use `var` or `std`')
    return {
        'res': u'[{0} - {1} - {2}'
        .format(conf_lower, mean, conf_upper),
        'u': conf_upper,
        'l': conf_lower,
        'ul': (conf_upper, conf_lower),
        'rng': conf_range,
        'mean': mean,
        'var': variance,
        'std': std,
        't': tstud
    }


def conf_clean(data):
    '''
    Removes data points that don't fit confidence level
    '''
    # TODO: add an option to choose between mask and data
    mask = np.array([], dtype=bool)
    for i in range(len(data)):
        data_miss = np.delete(data, i)
        (cu, cl) = conf_level(data_miss)['ul']
        if cl <= data[i] <= cu:
            mask = np.append(mask, [True])
        else:
            mask = np.append(mask, [False])
        return mask


def moving_average(array, size=50):
    '''
    return array with moving average points

    `size` amount of first points are `np.nan`

    Args:
        array: list type or something - input data
        size: int - amount of smoothing, points

    Returns:
        numpy.ndarray - smoothed data points
    '''
    mv_average = np.array([np.nan for _ in range(size)])
    for row in range(size, array.size):
        mv_average = np.append(mv_average, np.mean(array[row - size:row]))
    return mv_average


def sqdistance(ar1, ar2, order=1):
    '''
    calculates squared distance between points of two arrays
    '''
    try:
        iter(ar1)
        iter(ar2)
    except TypeError:
        raise TypeError("input arguments should be iterable!")
    if len(ar1) != len(ar2):  # Проверка на одинаковую длину
        raise TypeError("arguments should be same length!")
    dist = 0
    for item1, item2 in zip(ar1, ar2):
        dist += abs(item1**order - item2**order)
    return dist
