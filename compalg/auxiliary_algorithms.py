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


def poly_pth_root(f):
    """
    Computes the pth root of the polynomial f in F_q[x] where q = p^r, p prime and f' = 0 (not checked)
    There exists a polynomial g in F_q[x] such that f = g^p, i.e. a pth root of f, iff f' = 0.
    """

    field = f.parent()
    p = field.base().characteristic()
    q = field.base().order()

    coefs = f.coefficients(sparse=False)
    root_coefs = [0] * (1 + (len(coefs) / int(p)))

    for i, coef in enumerate(reversed(coefs)):
        if coef != 0:
            j = len(coefs) - i - 1
            root_coefs[j / p] = field(repeated_square(coef, q / p, q))  # a_j^{p^(r-1)} = a_j^{q/p}

    return field(root_coefs)


def poly_content(f):
    """
    Computes the polynomial content of f in Z[x]. The polynomial content is the gcd of its coefficients.
    """
    return f.parent().base()(gcd(map(ZZ, f.coefficients())))
