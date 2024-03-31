import numpy as np
from math import factorial
from scipy.stats import norm
from decimal import *
from termcolor import colored


#function gets name for each item in result
_add_names_for_result = lambda _inp: {'Bernoulli ' : _inp[0], 'Poisson ' : _inp[1], 'Local Moivre-Laplace ' : _inp[2], 'Integral Moivre-Laplace ' : _inp[3]}


# binomial coefficients counter
bin_coef = lambda n, k: factorial(n) / (factorial(k) * factorial(n - k))


# Bernoulli function
def bernoulli(k, p, n):
    return bin_coef(n, k) * p ** k * (1 - p) ** (n - k)


# Poisson function
def poisson(k, p, n):
    _lambda = p * n
    try:
      return Decimal(_lambda ** k * np.exp(-_lambda)) / Decimal(factorial(k))

    except OverflowError:
      print(colored('Attention! For this function probability can not be calculated by Poisson theorem!', 'red'))
      return "Can't count Poisson function!"


# integral Moivre--Laplace function
def integral_ML(k1, k2, p, n):
    return norm.cdf((k2 - n * p) / np.sqrt(n * p * (1 - p))) - norm.cdf((k1 - n * p) / np.sqrt(n * p * (1 - p)))


# local Moivre--Laplace function
def local_ML(k, p, n):
    return norm.pdf((k - n * p) / np.sqrt(n * p * (1 - p))) / np.sqrt(n * p * (1 - p))


# k in some area of n / 2
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
    _bern += bernoulli(k, p, n)
    _locM += local_ML(k, p, n)

  _intM = integral_ML(k1, k2, p, n)

  res = [_bern, _pois, _locM, _intM]

  return _add_names_for_result(res)


# Sn leq 5 fun
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
    _bern += bernoulli(k, p, n)
    _locM += local_ML(k, p, n)

  _intM = integral_ML(0, 5, p, n)

  res = [_bern, _pois, _locM, _intM]

  return _add_names_for_result(res)


# for k max
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
    _bern += bernoulli(k, p, n)
    _locM += local_ML(k, p, n)

  _intM = integral_ML(0, max_k, p, n)

  res = [_bern, _pois, _locM, _intM]

  return _add_names_for_result(res)


def tester():
    n = 1000
    p = 0.01

    leq_5 = k_leq_5_counter(n, p)

    print('For P(Sn leq 5) \n')

    for key, value in leq_5.items():
      print("{0}: {1}".format(key,value))

    print('\n ------------------- \n')

    print('For P(in area of n / 2) \n')

    in_area = k_in_area_counter(n, p)

    for key, value in in_area.items():
      print("{0}: {1}".format(key,value))

    print('\n ------------------- \n')

    print('For P(k*) \n')

    k_max = k_max_prob(n, p)

    for key, value in k_max.items():
      print("{0}: {1}".format(key,value))


tester()

