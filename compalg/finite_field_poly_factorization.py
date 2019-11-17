"""
Factoring algorithm of a polynomial in a finite field.
"""
from sage.all import *

from auxiliary_algorithms import repeated_square
from gcd_ufd import gcd_ufd


def squarefree_decomposition(f):
    """
    Square-free factorization.

    Computes the square-free decomposition of a nonconstant polynomial f in a finite field F_q[x]. A polynomial is
    called square-free if it is not divisible by the square of any polynomial of degree greater than zero.
    This decomposition is a sequence (f_1, ..., f_k) of polynomials in F_q[x] such that f = (f_1)^(s_1)*...*(f_k)^(s_k)

    Parameters
    ----------
    f : A nonconstant, monic polynomial in F_q[x].

    Returns
    -------
    The square-free decomposition [(f_1, s_1), ..., (f_k, s_k)] of f.
    """

    field = f.parent()
    base_field = field.base()

    if not base_field.is_field() or not base_field.is_finite() or not base_field.characteristic().is_prime():
        raise ValueError("The base field must be a finite field with a prime number of elements")

    if f.degree() < 1:
        return ValueError("f must be a nonconstant polynomial")
    if not f.is_monic():
        return ValueError("f must be a monic polynomial")

    p = base_field.characteristic()
    fs = []
    s = 1

    while f != field.one():
        j = 1
        f1 = f.diff()
        gcd_aux = gcd_ufd(f, f1)
        aux = f / gcd_aux
        g = field(f / gcd_aux)
        while g != field.one():
            f = field(f / g)
            h = gcd_ufd(f, g)
            m = field(g / h)
            if m != field.one():
                fs.append((m, j*s))
            g = h
            j += 1
        if f != field.one():  # f is a pth power
            f.pth_root()
            s = p*s

    return fs


def distinct_degree_decomposition(f):
    """
    Distinct-degree factorization.

    Computes the distinct-degree decomposition of a nonconstant polynomial f in a finite field Fq[x].
    This decomposition is a sequence (g_1, ..., g_s) of polynomials, where g_i is the product of all monic irreducible
    polynomials in F_q[x] of degree i that divide f, and f_s != 1.

    Parameters
    ----------
    f : A nonconstant, squarefree, monic polynomial in F_q[x].

    Returns
    -------
    The distinct-degree decomposition (g_1, ..., g_s) of f.
    """

    field = f.parent()
    base_field = field.base()

    if not base_field.is_field() or not base_field.is_finite() or not base_field.characteristic().is_prime():
        raise ValueError("The base field must be a finite field with a prime number of elements")

    if f.degree() < 1:
        return ValueError("f must be a nonconstant polynomial")
    if not f.is_monic():
        return ValueError("f must be a monic polynomial")
    if not f.is_squarefree():
        return ValueError("f must be a squarefree polynomial")

    q = base_field.order()

    h = field('x')
    f_i = f
    i = 0
    g = []

    while f_i.degree() >= 2*(i + 1):  # Early abort
        i += 1
        h = repeated_square(h, q, f)
        g_i = gcd_ufd(h - field('x'), f_i)
        g.append(g_i)
        f_i = field(f_i / g_i)

    if f_i != field.one():
        g.append(f_i)

    return g


def main():
    """ Execute the examples. """

    # Distinct degree decomposition
    F3 = GF(3)
    R = PolynomialRing(F3, 'x')
    f1 = R('x*(x + 1)*(x^2 + 1)*(x^2 + x + 2)')
    print distinct_degree_decomposition(f1)  # Expected: [x^2 + x, x^4 + x^3 + x + 2]

    f2 = R('x^8 + x^7 - x^6 + x^5 - x^3 - x^2 - x')
    print distinct_degree_decomposition(f2)  # Expected: [x, x^4 + x^3 + x - 1, x^3 - x + 1]

    f3 = R('(x^5 + x^2 + x + 1)^2 * x^5')  # Expected: [(x^5 + x^2 + x + 1, 2), (x, 5)]
    print squarefree_decomposition(f3)

    f4 = R('x^11 + 2*x^9 + 2*x^8 + x^6 + x^5 + 2*x^3 + 2*x^2 + 1')  # TODO: FIX
    print squarefree_decomposition(f4)


if __name__ == '__main__':
    main()
