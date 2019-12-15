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

    order = base_field.order()
    b = base_field.one()
    i = 0
    while b != alpha and i < order:
        b = gen * b
        i += 1

    if i == order:
        raise ValueError("No discrete log of %s found to base %s using gen %s" % (alpha, order, gen))

    return i


def discrete_log_bsgs(gen, alpha):
    """
    Baby step / giant step discrete logarithm computing.

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

    order = base_field.order()
    m = isqrt(order) + 1
    table = dict()

    # Baby steps, fill the table
    b = base_field.one()
    for i in xrange(m):
        table[b] = i
        b = b * gen

    # Giant steps
    genp = gen ** -m
    b = alpha
    for i in xrange(m):
        j = table.get(b)
        if j is not None:
            return i*m + j
        b = b * genp

    raise ValueError("No discrete log of %s found to base %s using gen %s" % (alpha, order, gen))


def main():
    K = GF(137)
    print discrete_log_bsgs(K(7), K(11))  # Expected 45

    K = GF(3 ** 6, 'b')
    b = K.gen()
    a = b ** 210
    print discrete_log_bsgs(b, a)  # Expected 210

    try:
        K = GF(37)
        print discrete_log_bsgs(K(1), K(2))  # No discrete log of 2 found to base 37 using gen 1
    except ValueError as err:
        print err


if __name__ == '__main__':
    main()
