"""
Factoring algorithm of a polynomial in a finite field.
"""
from sage.all import *

from auxiliary_algorithms import repeated_square, poly_pth_root
from gcd_ufd import gcd_ufd


def squarefree_decomposition(f):
    """
    Square-free factorization.

    Computes the square-free decomposition of a nonconstant monic polynomial f in a finite field F_q[x]. A polynomial is
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
        raise ValueError("The base field must be a finite field with a prime power number of elements")

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
        g = field(f / gcd_ufd(f, f1).monic())
        while g != field.one():
            f = field(f / g)
            h = gcd_ufd(f, g).monic()
            m = field(g / h)
            if m != field.one():
                fs.append((m, j * s))
            g = h
            j += 1
        if f != field.one():  # f is a pth power
            f = poly_pth_root(f)
            s = p * s

    return fs


def distinct_degree_decomposition(f):
    """
    Distinct-degree factorization.

    Computes the distinct-degree decomposition of a nonconstant, squarefree, monic polynomial f in a finite field
    F_q[x]. This decomposition is a sequence (g_1, ..., g_k) of polynomials, where g_i is the product of all monic
    irreducible polynomials in F_q[x] of degree s_i that divide f, and f_k != 1.

    Parameters
    ----------
    f : A nonconstant, squarefree, monic polynomial in F_q[x].

    Returns
    -------
    The distinct-degree decomposition ((g_1, s_1), ..., (g_k, s_k)) of f.
    """

    field = f.parent()
    base_field = field.base()

    if not base_field.is_field() or not base_field.is_finite() or not base_field.characteristic().is_prime():
        raise ValueError("The base field must be a finite field with a prime power number of elements")

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
        g_i = gcd_ufd(h - field('x'), f_i).monic()
        if g_i != field.one():
            g.append((g_i, i))
        f_i = field(f_i / g_i)

    if f_i != field.one():
        g.append((f_i, i + 1))

    return g


def __equal_degree_splitting(f, d):
    """
    Outputs a proper monic factor of f or None if it fails.
    Fails with probability less than 2^(1-r) where r = n/d >= 2.
    """

    field = f.parent()
    q = field.base().order()

    a = field.random_element()
    while a.degree() >= f.degree() or a.degree() < 1:
        a = field.random_element()

    g1 = gcd_ufd(a, f).monic()
    if g1 != field.one():
        return g1

    b = repeated_square(a, (q**d - 1) / 2, f)
    g2 = gcd_ufd(b - 1, f).monic()
    if g2 != field.one() and g2 != f:
        return g2
    else:
        return None


def __recursive_equal_degree_decomposition(f, d):
    """
    Recursive implementation of equal_degree_decomposition. Doesn't check input parameters.
    """

    field = f.parent()

    if f.degree() == d:
        return [f]

    g = __equal_degree_splitting(f, d)
    while g is None:
        g = __equal_degree_splitting(f, d)

    return __recursive_equal_degree_decomposition(g, d) + __recursive_equal_degree_decomposition(field(f / g), d)


def equal_degree_decomposition(f, d):
    """
    Equal-degree factorization.

    Computes the monic irreducible factors of of a squarefree, monic polynomial f of degree n > 0 in a finite field
    F_q[x], where q is an odd prime power and a divisor d of n, so that all irreducible factors of f have degree d.

    Parameters
    ----------
    f : A quarefree, monic polynomial of degree n > 0 in F_q[x], where q is an odd prime power.
    d: A divisor of n, where n is the degree of f, so that all irreducible factors of f have degree d.

    Returns
    -------
    The monic irreducible factors of f in F_q[x].
    """

    field = f.parent()
    base_field = field.base()

    if not base_field.is_field() or not base_field.is_finite() or not base_field.characteristic().is_prime() \
            or base_field.characteristic() == 2:
        raise ValueError("The base field must be a finite field with an odd prime power number of elements")

    if f.degree() < 1:
        return ValueError("f must be a nonconstant polynomial")
    if not f.is_monic():
        return ValueError("f must be a monic polynomial")
    if not f.is_squarefree():
        return ValueError("f must be a squarefree polynomial")

    if f.degree() % d != 0:
        return ValueError("d must be a divisor of degree(f)")

    return __recursive_equal_degree_decomposition(f, d)


def poly_factorization(f):
    """
    Polynomial factorization in a finite field

    Factors a nonconstant monic polynomial f in a finite field F_q[x], where q is an odd prime power q = p^r.

    Parameters
    ----------
    f : A nonconstant, monic polynomial in F_q[x].

    Returns
    -------
    The factorization [(f_1, s_1), ..., (f_k, s_k)] of f.
    """

    irreducible_factors = []
    for squarefree, power in squarefree_decomposition(f):
        for distinct_degree, d in distinct_degree_decomposition(squarefree):
            for irreducible_factor in equal_degree_decomposition(distinct_degree, d):
                irreducible_factors.append((irreducible_factor, power))

    return irreducible_factors


def main():
    """ Execute the examples. """

    F3 = GF(3)
    R = PolynomialRing(F3, 'x')

    f1 = R('(x^5 + x^2 + x + 1)^2 * x^5')
    print squarefree_decomposition(f1)  # Expected: [(x^5 + x^2 + x + 1, 2), (x, 5)]

    f2 = R('x^11 + 2*x^9 + 2*x^8 + x^6 + x^5 + 2*x^3 + 2*x^2 + 1')
    print squarefree_decomposition(f2)  # Expected: [(x + 1, 1), (x + 2, 4), (x^2 + 1, 3)]

    f3 = R('x*(x + 1)*(x^2 + 1)*(x^2 + x + 2)')
    print distinct_degree_decomposition(f3)  # Expected: [(x^2 + x, 1), (x^4 + x^3 + x + 2, 2)]

    f4 = R('x^8 + x^7 - x^6 + x^5 - x^3 - x^2 - x')
    print distinct_degree_decomposition(f4)  # Expected: [(x, 1), (x^4 + x^3 + x + 2, 2), (x^3 + 2*x + 1, 3)]

    f5 = R('x^4 + x^3 + x - 1')
    print equal_degree_decomposition(f5, 2)  # Expected: [x^2 + 1, x^2 + x + 2]

    f6 = R('x^9 + x^8 - x^7 + x^6 - x^4 - x^3 - x^2')
    print poly_factorization(f6)  # Expected: [(x^2 + x + 2, 1), (x^2 + 1, 1), (x^3 + 2*x + 1, 1), (x, 2)]


if __name__ == '__main__':
    main()
