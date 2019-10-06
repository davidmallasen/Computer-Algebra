"""
Extended euclidean algorithm (gcd).
"""
import operator
from sage.all import *


def extended_euclidean_algorithm(f, g):
    """
    Euclidean Algorithm.

    Returns one gcd of two elements f, g belonging to an euclidean domain.

    Parameters
    ----------
    f : the first element
    g : the second element

    Returns
    -------
    r : a list whose nth element is the remainder of the nth step of the algorithm.
    s : a list that verifies that s[i] * f + t[i] * g = r[i].
    t : a list that verifies that s[i] * f + t[i] * g = r[i].
    q : a list whose nth element is the quotient of the nth step of the algorithm. q[0] is a dummy element.
    """

    if f.parent() is not g.parent():
        raise ValueError("Arguments should belong to the same ring")

    domain = f.parent()

    if not domain.is_euclidean_domain():
        raise ValueError("Arguments should belong to an euclidean domain")

    r = [f, g]
    s = [domain.one(), domain.zero()]
    t = [domain.zero(), domain.one()]
    q = [domain.zero()]

    i = 1
    while r[i] != domain.zero():
        quo, _ = r[i - 1].quo_rem(r[i])
        q.append(quo)
        r.append(r[i - 1] - q[i] * r[i])
        s.append(s[i - 1] - q[i] * s[i])
        t.append(t[i - 1] - q[i] * t[i])
        i += 1

    return r, s, t, q


def result_to_pretty_string(f, g, rs, ss, ts):
    """
    Creates a nicely formatted string using the result from the algorithm.
    """
    return '\n'.join(map(lambda (r, s, t): '{r} = {s} * {f} + {t} * {g}'.format(r=r, s=s, t=t, f=f, g=g), zip(rs, ss, ts)))


def main():
    """ Execute the examples. """
    print extended_euclidean_algorithm(ZZ(126), ZZ(35))

    R = PolynomialRing(QQ, 'x')
    print extended_euclidean_algorithm(R('x^3 - 1'), R('x + 1'))


if __name__ == '__main__':
    main()
