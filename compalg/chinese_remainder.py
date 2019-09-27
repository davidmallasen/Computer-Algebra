"""
Chinese remainder theorem algorithm (compute the inverse).
"""
import operator
from sage.all import ZZ


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

    parent_ = residues[0].parent()
    try:
        if not parent_.is_euclidean_domain():
            raise TypeError("All the elements of both lists must be in an Euclidean Domain")
    except AttributeError:  # In case parent_ doesnt have the is_euclidean_domain method
        raise TypeError("All the elements of both lists must be in an Euclidean Domain")

    type_ = type(residues[0])
    for i, _ in enumerate(residues):
        if not isinstance(residues[i], type_) or not isinstance(moduli[i], type_):
            raise TypeError("All the elements of both lists must be in the same Euclidean Domain")

    # Compute the actual solution
    m = reduce(operator.mul, moduli)
    c = 0
    for (v_i, m_i) in zip(residues, moduli):
        quo, _ = m.quo_rem(m_i)
        # ..., s_i, ... = extended_euclid(..., quo, ...)    TODO
        s_i = 0  # TODO
        _, c_i = (v_i * s_i).quo_rem(m_i)
        c += c_i * quo
    return c


def main():
    """ Execute the examples. """
    print chinese_remainder([ZZ(2), ZZ(7)], [ZZ(11), ZZ(13)])  # Sol = 46


if __name__ == '__main__':
    main()
