"""
Euclidean algorithm (gcd).
"""
import operator
from sage.all import *


def euclidean_algorithm(f, g):
    """
    Euclidean Algorithm.

    Returns one gcd of two elements f, g belonging to an euclidean domain.

    Parameters
    ----------
    f : the first element
    g : the second element

    Returns
    -------
    A gcd of f and g.
    """

    if f.parent() is not g.parent():
        raise ValueError("Arguments should belong to the same ring")

    domain = f.parent()

    if not domain.is_euclidean_domain():
        raise ValueError("Arguments should belong to an euclidean domain")

    a = f
    b = g

    while b != domain.zero():
        _, rem = a.quo_rem(b)
        a, b = b, rem

    return a


def main():
    """ Execute the examples. """
    print euclidean_algorithm(ZZ(126), ZZ(35))  # Expected 7
    print euclidean_algorithm(ZZ(112242), ZZ(989712))  # Expected 6
    print euclidean_algorithm(ZZ(33461), ZZ(236322))  # Expected 1

    R = PolynomialRing(QQ, 'x')
    print euclidean_algorithm(R('x^2 - 1'), R('x + 1'))  # Expected x+1
    print euclidean_algorithm(R('x + 1'), R('x + 2'))  # Expected 1 (or -1)


if __name__ == '__main__':
    main()
