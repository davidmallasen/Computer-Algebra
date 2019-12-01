"""
Ideal membership algorithm. Does an element belong to a given ideal?
"""
from sage.all import *
from buchberger_algorithm import buchberger_algorithm
from multivariate_division import multivariate_division_with_remainder
from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
from sage.rings.polynomial.multi_polynomial_ring import is_MPolynomialRing
from sage.rings.ideal import is_Ideal


def is_member(f, I):
    """
    Determines if the polynomial f is a member of the ideal I.

    Parameters
    ----------
    f: a polynomial of a multivariate polynomial ring over a field.
    I: an ideal of the same polynomial ring.

    Returns
    -------
    True if f is a member of the ideal I, False otherwise.
    """
    if not is_Ideal(I):
        raise TypeError('Argument should be an ideal')

    poly_ring = I.ring()

    if not f.parent() == poly_ring:
        raise ValueError('f must belong to the same polynomial ring that the ideal belongs to')

    if (not is_PolynomialRing(poly_ring)) and (not is_MPolynomialRing(poly_ring)):
        raise TypeError('The ideal should be of a polynomial ring')

    base_field = I.base_ring()

    if not base_field.is_field():
        raise TypeError('The ideal should be of a polynomial ring over a field')

    # Computation of a groebner basis
    G = buchberger_algorithm(I)

    # Since groebner basis have a unique remainder with multivariate division, f belongs to I <=> f rem G == 0
    # Formalized on theorem 21.28 of Modern Computer Algebra
    _, rem = multivariate_division_with_remainder(f, G)
    return rem == 0


def main():
    """ Execute the examples. """
    # Pag 609 modern computer algebra
    poly_ring = PolynomialRing(QQ, 'x,y,z', order='deglex')     # You can specify, e.g., order='lex' to override <
    x, y = poly_ring('x'), poly_ring('y')
    f1 = x ** 3 - 2 * x * y
    f2 = x ** 2 * y - 2 * y ** 2 + x
    I = poly_ring.ideal(f1, f2)

    # Examples generated with sage 'f in I'
    print is_member(f1, I)      # Expected true
    print is_member(f2, I)      # Expected true
    print is_member(f2 * x ** 3, I)     # Expected true
    print is_member(f1 + f2, I)         # Expected true
    print is_member(f1 * f2 ** 2, I)    # Expected true
    print is_member(f1 + f2, I)         # Expected true

    print is_member(x, I)    # Expected false
    print is_member(y, I)    # Expected false
    print is_member(poly_ring(1), I)    # Expected false
    print is_member(x + y, I)    # Expected false


if __name__ == '__main__':
    main()