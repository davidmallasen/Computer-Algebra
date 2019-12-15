"""
Multivariate division with remainder.
"""
from sage.all import *
from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
from sage.rings.polynomial.multi_polynomial_ring import is_MPolynomialRing


def multivariate_division_with_remainder(f, fs):
    """
    Performs multivariate division with remainder, that is, returns quotients q_1, ..., q_s and a remainder r such that
    q_1 * f_1 + ... + q_s * f_s + r = f, and no monomial in r is divisible by any of lt(f_i). Both f and all fs must
    belong to the same multivariate polynomial ring.

    Parameters
    ----------
    f : the numerator.
    fs : the list of denominators.

    Returns
    -------
    qs: the quotients.
    r: the remainder, f rem (f1, ..., fs).
    """
    poly_ring = f.parent()

    if (not is_PolynomialRing(poly_ring)) and (not is_MPolynomialRing(poly_ring)):
        raise TypeError('f and the fs should belong to a polynomial ring')

    base_field = poly_ring.base_ring()

    if not base_field.is_field():
        raise TypeError('f and the fs should belong to a polynomial ring over a field')

    if not all([poly_ring == g.parent()] for g in fs):
        raise ValueError("All polynomials must belong to the same ring")

    r = poly_ring(0)
    p = f
    q = [poly_ring(0)] * len(fs)

    while not p.is_zero():
        any_divides = False

        for i in range(len(fs)):
            fi = fs[i]

            if fi.lt().divides(p.lt()):
                div, _ = p.lt().quo_rem(fi.lt())
                q[i] += div
                p -= div * fi

                any_divides = True
                break

        if not any_divides:
            r += p.lt()
            p -= p.lt()

    return q, r


def main():
    """ Execute the examples. """
    poly_ring = PolynomialRing(QQ, 'x,y', order='lex')
    x, y = poly_ring('x'), poly_ring('y')
    f = x * y ** 2 - x
    f1 = x * y + 1
    f2 = y ** 2 - 1
    print multivariate_division_with_remainder(f, [f1, f2])     # Expected [y, 0], -x-y


if __name__ == '__main__':
    main()
