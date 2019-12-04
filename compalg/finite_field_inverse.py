"""
Computation of the inverse in a finite field.
"""
from extended_euclidean_algorithm import normalized_extended_euclidean_algorithm
from sage.all import *


def inverse_element(a):
    """
    Inverse in a finite field.

    Parameters
    ----------
    a : the element to be inverted, belonging to a Galois field.

    Returns
    -------
    Its inverse.
    """
    ring = a.parent()
    polynomial = ring.polynomial()

    return ring(__inverse_element(a.polynomial(), polynomial))


def __inverse_element(a, f_mod):
    """
    Inverse in a finite field.

    Returns the inverse of the given element.

    Parameters
    ----------
    a : the element to be inverted.
    f_mod : the polynomial to calculate the inverse with. We will calculate a^-1 in base_field / f_mod. f_mod must be
    irreducible.

    Both arguments must belong to the same finite field, with a prime number of elements.

    Returns
    -------
    Its inverse
    """

    base_field = a.parent().base()

    if not a.parent() == f_mod.parent():
        raise ValueError("Both polynomials must belong to the base_field's polynomial field")

    if not base_field.is_field() or not base_field.is_finite() or not base_field.characteristic().is_prime():
        raise ValueError("The base field must be a finite field with a prime number of elements")

    _, _, t, _ = normalized_extended_euclidean_algorithm(f_mod, a)
    return t[-2]


def main():
    """ Execute the examples. """
    Z5 = GF(5)
    R = PolynomialRing(Z5, 'x')
    f = R('x^3 - x + 2')
    a = R('x^2')
    print(__inverse_element(a, f))    # Expected x^2 + 2x - 1

    Z11 = GF(11)
    R = PolynomialRing(Z11, 'x')
    f = R('x')
    a = R(5)
    print __inverse_element(a, f)

    Z9 = GF(9)
    a = Z9(Z9.variable_name())
    expected = Z9.one().quo_rem(a)[0]
    print 'Expected {e}, got {a}'.format(e=expected, a=inverse_element(a))


if __name__ == '__main__':
    main()
