"""
Chinese remainder theorem algorithm (compute the inverse).
"""
import operator
from sage.all import *

from extended_euclidean_algorithm import extended_euclidean_algorithm


def chinese_remainder(residues, moduli):
    """
    Chinese Remainder Algorithm (CRA).

    Return a solution to a Chinese Remainder Theorem problem. Given an Euclidean Domain ED and
    two lists of elements in ED of the same length, ''residues'' and ''moduli'', returns an element
    of DE such that this element is congruent with residues[i] mod moduli[i] for all i.

    Parameters
    ----------
    residues : list
        List of residues as stated above.

    moduli : list
        List of pairwise coprime moduli as stated above.

    Returns
    -------
    A solution to the Chinese Remainder Theorem for ''residues'' and ''moduli''.
    If the lists are empty, returns ZZ(0).

    """
    # Check input parameters, if the lists are empty, return ZZ(0)
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

    # Compute the actual solution
    m = reduce(operator.mul, moduli)
    c = domain.zero()
    for (v_i, m_i) in zip(residues, moduli):
        quo, _ = m.quo_rem(m_i)
        _, s, _, _ = extended_euclidean_algorithm(quo, m_i)
        s_i = s[-2]     # s_i * quo + t_i * m_i = 1
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
    _example([R('2'), R('3')], [f, g])  # Not working!
    print crt([R('2'), R('3')], [f, g])  # Sage prebuilt function


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
