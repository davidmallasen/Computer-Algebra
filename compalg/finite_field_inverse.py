"""
Computation of the inverse in a finite field.
"""
from extended_euclidean_algorithm import extended_euclidean_algorithm
from sage.all import *


def inverse_element(a, f_mod):
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

    r, _, t, _ = extended_euclidean_algorithm(f_mod, a)
    return t[-2] * r[-2]   # Product by r[-2] to get the inverse if we calculate gcd() = -1


def main():
    """ Execute the examples. """
    Z5 = GF(5)
    R = PolynomialRing(Z5, 'x')
    f = R('x^3 - x + 2')
    a = R('x^2')
    print(inverse_element(a, f))    # Expected x^2 + 2x - 1


if __name__ == '__main__':
    main()
