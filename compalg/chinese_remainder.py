"""
Chinese remainder theorem algorithm (compute the inverse).
"""
from sage.all import *

from extended_euclidean_algorithm import normalized_extended_euclidean_algorithm


def positive_bezout_coefficients(a, b, normal=None):
    """
    Return x, y such that x * a + y * b = 1, assuming that a and b are coprime. Will throw an error if a and b are not
    coprime.

    Parameters
    ----------
    a : the first element.
    b : the second element.
    normal : a function R -> R, that returns a normal form for a given element (see normalized extended euclidean
            algorithm for more info).

    Returns
    -------
    x, y such that x * a + y * b = 1.
    """

    domain = a.parent()
    r, s, t, _ = normalized_extended_euclidean_algorithm(a, b, normal)
    gcd = r[-2]
    if gcd == domain.one():
        return s[-2], t[-2]
    else:
        raise ValueError('Arguments should be coprime, but their gcd is {gcd}'.format(gcd=gcd))


def chinese_remainder(residues, moduli, normal=None):
    """
    Chinese Remainder Algorithm (CRA).

    Return a solution to a Chinese Remainder Theorem problem. Given an Euclidean Domain ED and two lists of elements in
    ED of the same length, ''residues'' and ''moduli'', returns an element in ED such that this element is congruent
    with residues[i] mod moduli[i] for all i.

    Parameters
    ----------
    residues : a list of residues as stated above.
    moduli : a list of pairwise coprime moduli as stated above.
    normal : a function R -> R, that returns a normal form for a given element (see normalized extended euclidean
            algorithm for more info).

    Returns
    -------
    A solution to the Chinese Remainder Theorem for ''residues'' and ''moduli''.
    If the lists are empty, returns ZZ(0).
    """

    if not isinstance(residues, list) or not isinstance(moduli, list):
        raise TypeError("Arguments should be lists")
    if len(residues) != len(moduli):
        raise ValueError("Arguments should be lists of the same length")

    if not residues:
        return ZZ(0)

    domain = residues[0].parent()
    try:
        if not domain.is_euclidean_domain():
            raise TypeError("All the elements of both lists must be in an Euclidean Domain")
    except AttributeError:  # In case parent_ doesnt have the is_euclidean_domain method
        raise TypeError("All the elements of both lists must be in an Euclidean Domain")

    type_ = type(residues[0])
    for i, _ in enumerate(residues):
        if not isinstance(residues[i], type_) or not isinstance(moduli[i], type_):
            raise TypeError("All the elements of both lists must be in the same Euclidean Domain")

    m = prod(moduli)
    c = domain.zero()
    for (v_i, m_i) in zip(residues, moduli):
        quo, _ = m.quo_rem(m_i)
        s_i, _ = positive_bezout_coefficients(quo, m_i, normal)   # s_i * quo + t_i * m_i = 1
        _, c_i = (v_i * s_i).quo_rem(m_i)
        c += c_i * quo
    return c


def main():
    """ Execute the examples. """
    _example([ZZ(2), ZZ(7)], [ZZ(11), ZZ(13)])  # Expected 46
    _example([ZZ(2), ZZ(3), ZZ(2)], [ZZ(3), ZZ(5), ZZ(7)])  # Expected 128

    R = PolynomialRing(QQ, 'x')
    f = R('x - 1')
    g = R('x - 2')
    _example([R('2'), R('3')], [f, g])  # Expected x + 1


def _example(residues, moduli):
    """ Executes and verifies an example. Prints the output."""
    ans = chinese_remainder(residues, moduli)
    correct = true
    for r, m in zip(residues, moduli):
        _, rem = ans.quo_rem(m)
        if rem != r:
            correct = false
            break
    print correct, ans


if __name__ == '__main__':
    main()
