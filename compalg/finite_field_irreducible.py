"""
Irreducibility check in a finite field.
"""
from extended_euclidean_algorithm import extended_euclidean_algorithm
from sage.all import *


def __repeated_square(a, n):
    """
    Computes a^n, with a belonging to a ring with 1 (not checked) and n > 0.
    """
    ns = list(bin(n)[2:])[1:]

    b = a
    for coef in ns:
        print(b)
        if coef == '1':
            b = b * b * a
        else:
            b = b * b

    return b

def inverse_element(f):
    """
    Irreducibility check for polynomials in a finite field.

    Parameters
    ----------
    f : the element to be checked. Must be a polynomial with coefficients belonging to a finite field.

    Returns
    -------
    A boolean representing whether f is irreducible or not.
    """

    poly_field = f.parent()
    base_field = poly_field.base()

    if not base_field.is_field() or not base_field.is_finite():
        raise ValueError("The base field must be a finite field")

    n = f.degree()
    x_power = __repeated_square(R(R.variable_name()), n)

    # Algorithm in page 423 of MCA


def main():
    """ Execute the examples. """
    Z5 = GF(5)
    R = PolynomialRing(Z5, 'x')
    f = R('x^3 - x + 2')
    a = R('x^2')
    print __repeated_square(ZZ(5), 3)


if __name__ == '__main__':
    main()
