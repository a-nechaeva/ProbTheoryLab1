import numpy as np
from math import factorial
from scipy.stats import norm
from decimal import *
from termcolor import colored

#  function gets name for each item in result
_add_names_for_result = lambda _inp: {'Bernoulli ': _inp[0], 'Poisson ': _inp[1], 'Local Moivre-Laplace ': _inp[2],
                                      'Integral Moivre-Laplace ': _inp[3]}


#  for smart print
def _smart_print(res):
    for key, value in res.items():
        print('\u001b[33;1m', key, '\u001b[0m \n')
        for kkey, vvalue in value.items():
            if isinstance(vvalue, str):
                print('\u001b[31;1m Attention! Probability can not be calculated by ', kkey, '! \u001b[0m ')
            else:
                print("\u001b[34;1m {0} \u001b[0m : \u001b[32;1m {1} \u001b[0m".format(kkey, vvalue))
        print('\u001b[36m ------------------------------------------------------------------------------------- \n')


#  binomial coefficients counter
bin_coef = lambda n, k: Decimal(factorial(n)) / Decimal((factorial(k) * factorial(n - k)))


#  Bernoulli function
def bernoulli(k, p, n):
    try:
        return bin_coef(n, k) * Decimal(p ** k * (1 - p) ** (n - k))
    except OverflowError:
        return "Can't count Bernoulli function!"


#  Poisson function
def poisson(k, p, n):
    _lambda = p * n
    try:
        return Decimal(_lambda ** k * np.exp(-_lambda)) / Decimal(factorial(k))

    except OverflowError:
        return "Can't count Poisson function!"


#  integral Moivre--Laplace function
def integral_ML(k1, k2, p, n):
    return norm.cdf((k2 - n * p) / np.sqrt(n * p * (1 - p))) - norm.cdf((k1 - n * p) / np.sqrt(n * p * (1 - p)))


#  local Moivre--Laplace function
def local_ML(k, p, n):
    return norm.pdf((k - n * p) / np.sqrt(n * p * (1 - p))) / np.sqrt(n * p * (1 - p))


#  k in some area of n / 2
def k_in_area_counter(n, p):
    magic_sqrt = np.sqrt(n * p * (1 - p))
    k1 = int(n / 2 - magic_sqrt)
    k2 = int(n / 2 + magic_sqrt)

    _bern = 0
    _pois = 0
    _locM = 0

    for k in range(k1, k2 + 1):
        t = poisson(k, p, n)

        if isinstance(t, str):
            _pois = t
            break

        _pois += t

    for k in range(k1, k2 + 1):
        t = bernoulli(k, p, n)

        if isinstance(t, str):
            _bern = t
            break

        _bern += t

    for k in range(k1, k2 + 1):
        _locM += local_ML(k, p, n)

    _intM = integral_ML(k1, k2, p, n)

    res = [_bern, _pois, _locM, _intM]

    return _add_names_for_result(res)


#  Sn leq 5 fun
def k_leq_5_counter(n, p):
    _bern = 0
    _pois = 0
    _locM = 0

    for k in range(6):
        t = poisson(k, p, n)

        if isinstance(t, str):
            _pois = t
            break

        _pois += t

    for k in range(6):
        t = bernoulli(k, p, n)

        if isinstance(t, str):
            _bern = t
            break

        _bern += t

    for k in range(6):
        _locM += local_ML(k, p, n)

    _intM = integral_ML(0, 5, p, n)

    res = [_bern, _pois, _locM, _intM]

    return _add_names_for_result(res)


#  for k max
def k_max_prob(n, p):
    max_k = int(n * p)
    _bern = 0
    _pois = 0
    _locM = 0

    for k in range(max_k):
        t = poisson(k, p, n)

        if isinstance(t, str):
            _pois = t
            break

        _pois += t

    for k in range(max_k):
        t = bernoulli(k, p, n)

        if isinstance(t, str):
            _bern = t
            break

        _bern += t

    for k in range(max_k):
        _locM += local_ML(k, p, n)

    _intM = integral_ML(0, max_k, p, n)

    res = [_bern, _pois, _locM, _intM]

    return _add_names_for_result(res)


def tester():
    n = [100, 1000, 10000]
    p = [0.001, 0.01, 0.1, 0.25, 0.5]

    for i in range(len(n)):
        for j in range(len(p)):
            print('', i * len(p) + j + 1, ': n = ', n[i], ' p = ', p[j], '')
            leq_5 = k_leq_5_counter(n[i], p[j])

            in_area = k_in_area_counter(n[i], p[j])

            k_max = k_max_prob(n[i], p[j])

            res = {'For P(Sn leq 5)': leq_5, 'For P(in area of n / 2)': in_area, 'For P(k*)': k_max}
            _smart_print(res)


tester()

