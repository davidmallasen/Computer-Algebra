"""
Discrete logarithm algorithms
"""
from sage.all import *


def discrete_log_brute(gen, alpha):
    """
    Brute force discrete logarithm computing.

    Return a discrete logarithm of alpha using gen as a generator of the finite field.

    Parameters
    ----------
    gen : Generator of the finite field.
    alpha : Element of the field.

    Returns
    -------
    A discrete logarithm of alpha.
    """

    base_field = gen.parent()

    if not gen.parent() == alpha.parent():
        raise ValueError("Both elements must belong to the same finite field")

    if not base_field.is_field() or not base_field.is_finite() or not base_field.characteristic().is_prime():
        raise ValueError("The base field must be a finite field with a prime number of elements")

    b = base_field.one()
    i = 0
    while b != alpha:
        b = gen * b
        i += 1
    return i


def main():
    Z137 = GF(137)
    print discrete_log_brute(Z137(7), Z137(11))  # Expected 45


if __name__ == '__main__':
    main()
