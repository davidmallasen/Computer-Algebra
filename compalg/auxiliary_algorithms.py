"""
Some auxiliary algorithms.
"""
from sage.all import *


def repeated_square(a, n, mod=None):
    """
    Computes a^n % mod, with a belonging to a ring with 1 (not checked) and n > 0.
    """
    ns = list(bin(n)[2:])[1:]

    b = a
    for coef in ns:
        if coef == '1':
            b = (b * b * a)
        else:
            b = (b * b)

        if mod is not None:
            b %= mod

    return b


def horners_evaluation(f, u):
    """
    Fast evaluation of f(u).
    """
    if f.is_zero():
        return f.parent().base().zero()

    coefs = reversed(f.coefficients(sparse=False))
    return reduce(lambda acc, next: acc * u + next, coefs)