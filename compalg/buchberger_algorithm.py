"""
Buchberger algorithm.
"""
from multivariate_division import multivariate_division_with_remainder
from sage.all import *
from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
from sage.rings.polynomial.multi_polynomial_ring import is_MPolynomialRing
from sage.rings.ideal import is_Ideal
import operator


def __multidegree(f):
    """
    Computes the multidegree of f.
    """

    mon = f.monomials()
    mon.sort()      # Sorts them with the ring order
    return mon[-1].degrees()


def __s_polynomial(g, h):
    """
    Computes the S-polynomial of g, h.
    """

    deg_g = __multidegree(g)
    deg_h = __multidegree(h)
    max_deg = map(max, zip(deg_g, deg_h))
    R = g.parent()

    # Builds a polynomial with the variables raised to max_deg, in order
    vars = map(R, R.variable_names())
    x_pow_max_deg = reduce(operator.mul, [x ** d for (d, x) in zip(max_deg, vars)], R(1))

    quo_g, _ = x_pow_max_deg.quo_rem(g.lt())
    quo_h, _ = x_pow_max_deg.quo_rem(h.lt())
    return quo_g * g - quo_h * h


def __unordered_pairs(l):
    """
    Generates a list of pairs from the given list, excluding reorders (if (a, b) is returned, (b, a) isn't), and
    excluding pairs where the two elements are the same.

    If two elements in the list are equal, this method will still consider them as different.
    """

    return [(l[i], l[j]) for i in range(len(l) - 1) for j in range(i + 1, len(l))]


def buchberger_algorithm(I):
    """
    Buchberger algorithm for the computation of a Groebner basis.

    Returns a list of polynomials that form a Groebner basis of the given ideal.

    Parameters
    ----------
    I: an ideal of a multivariate polynomial ring over a field.

    Returns
    -------
    The list of polynomials forming a Groebner basis.
    """
    if not is_Ideal(I):
        raise TypeError('Argument should be an ideal')

    poly_ring = I.ring()

    if (not is_PolynomialRing(poly_ring)) and (not is_MPolynomialRing(poly_ring)):
        raise TypeError('The ideal should be of a polynomial ring')

    base_field = I.base_ring()

    if not base_field.is_field():
        raise TypeError('The ideal should be of a polynomial ring over a field')

    # Initialize G with the generators of the ideal
    G = list(I.basis)

    # It will finish by construction of the algorithm
    while True:
        S = []
        G.sort(reverse=True)    # Sorts with the polynomial ring ordering

        for (g1, g2) in __unordered_pairs(G):
            r = __s_polynomial(g1, g2)
            _, r = multivariate_division_with_remainder(r, G)

            if r != 0:
                S.append(r)

        # If S is empty, we already have a Groebner basis
        if not S:
            return G

        G.extend(S)


def __test_s_polynomial():
    """
    Pag. 606 modern computer algebra
    """
    poly_ring = PolynomialRing(QQ, 'x,y', order='deglex')
    x, y = poly_ring('x'), poly_ring('y')
    g = x ** 3 - 2 * x * y
    h = x ** 2 * y - 2 * y ** 2 + x
    print __s_polynomial(g, h)      # Expected -x^2


def main():
    """ Execute the examples. """
    # Pag 609 modern computer algebra
    poly_ring = PolynomialRing(QQ, 'x,y,z', order='deglex')     # You can specify, e.g., order='lex' to override <
    x, y = poly_ring('x'), poly_ring('y')
    f1 = x ** 3 - 2 * x * y
    f2 = x ** 2 * y - 2 * y ** 2 + x
    I = poly_ring.ideal(f1, f2)
    print buchberger_algorithm(I)   # Expected f1, f2, -x^2, -2xy, -2y^2+x


if __name__ == '__main__':
    main()
