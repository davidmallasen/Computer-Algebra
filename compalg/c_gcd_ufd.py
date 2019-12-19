"""
GCD computing for members of any UCD.
"""
from sage.all import *
from sage.rings.polynomial.polynomial_ring import is_PolynomialRing

from a1_euclidean_algorithm import euclidean_algorithm


def __poly_content(f):
    """
    Computes the polynomial content for the given polynomial.
    """

    return reduce(euclidean_algorithm, f.coefficients())


def __gcd_ufd_primitive_polynomial(f, g):
    """
    GCD computing.

    Returns one gcd of two elements f, g belonging to a ring whose base ring is an Euclidean domain for f, g primitive
    polynomials. Doesn't check for consistency, this should happen in the caller method.

    Parameters
    ----------
    f : the first element, must be primitive.
    g : the second element, must be primitive.

    Returns
    -------
    The gcd.
    """

    # We convert the polynomials to polynomials in the fraction field. Since polynomials over a field are an euclidean
    # domain, we can apply the euclidean algorithm to get the gcd of the monic polynomials and the gcd of the leading
    # coefficients. We then adjust the result to be primitive, since the gcd of primitive polynomials is primitive.

    domain = f.parent()
    base_domain = domain.base()
    fraction_field_base = base_domain.fraction_field()
    polys_over_fraction_field = PolynomialRing(fraction_field_base, 'x')

    f_in_ff, g_in_ff = polys_over_fraction_field(f), polys_over_fraction_field(g)
    v = euclidean_algorithm(f_in_ff, g_in_ff).monic()
    b = euclidean_algorithm(f.leading_coefficient(), g.leading_coefficient())

    # h should be now an element of the original domain again, we don't want it as element of the fraction field, since
    # that breaks the polynomial content calculation
    h = domain(b * v)
    return h / __poly_content(h)


def gcd_ufd(f, g):
    """
    GCD computing.

    Returns one gcd of two elements f, g belonging to a ring whose base ring is an Euclidean domain.

    Parameters
    ----------
    f : the first element.
    g : the second element.

    Returns
    -------
    The gcd.
    """

    if f.parent() is not g.parent():
        raise ValueError("Arguments should belong to the same ring")

    domain = f.parent()

    if not is_PolynomialRing(domain):
        raise ValueError("Arguments should be polynomials")

    base_domain = domain.base()

    if not base_domain.is_euclidean_domain():
        raise ValueError("The base ring for the polynomial ring must be an Euclidean domain")

    # gcd(f, 0) = f, gcd(0, g) = g
    if f == domain.zero():
        if g == domain.zero():
            return domain.zero()
        else:
            return g
    else:
        if g == domain.zero():
            return f

    # We use the well known equalities: c(gcd(f, g)) = gcd(c(f), c(g)), pp(gcd(f, g)) = gcd(pp(f), pp(g)).

    cont_f, cont_g = __poly_content(f), __poly_content(g)
    cont_result = euclidean_algorithm(cont_f, cont_g)

    pp_f, _ = f.quo_rem(cont_f)
    pp_g, _ = g.quo_rem(cont_g)
    pp_result = __gcd_ufd_primitive_polynomial(pp_f, pp_g)

    return cont_result * pp_result


def main():
    """ Execute the examples. """
    R = PolynomialRing(QQ, 'x')
    print gcd_ufd(R('x^3 - 1'), R('x - 1'))     # Expected x - 1

    S = PolynomialRing(PolynomialRing(GF(5), 'y'), 'x')
    f = S('(y^3 + 3 * y^2 + 2 * y) * x^3 + (y^2 + 3 * y + 2) * x^2 + (y^3 + 3 * y^2 + 2 * y) * x + (y^2 + 3 * y + 2)')
    g = S('(2 * y^3 + 3 * y^2 + y) * x^2 + (3 * y^2 + 4 * y + 1) * x + (y + 1)')
    print gcd_ufd(f, g)    # Expected (y^2 + y) * x + (y + 1)


if __name__ == '__main__':
    main()
