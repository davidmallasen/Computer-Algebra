"""
Extended euclidean algorithm (gcd).
"""
from sage.all import *
from sage.rings.polynomial.polynomial_ring import is_PolynomialRing


def extended_euclidean_algorithm(f, g):
    """
    Extended euclidean Algorithm.

    Parameters
    ----------
    f : the first element.
    g : the second element.

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


def __lu(a, normal):
    """
    The leading unit of a, that is, the unit u such as a = u * normal(a).
    """

    if a.is_zero():
        return a.parent().one()

    quo, _ = a.quo_rem(normal(a))
    return quo


def normalized_extended_euclidean_algorithm(f, g, normal=None):
    """
    Normalized extended euclidean Algorithm.

    Parameters
    ----------
    f : the first element, belonging to the euclidean domain R.
    g : the second element, belonging to the euclidean domain R.
    normal : a function R -> R, that returns a normal form for a given element. That is, given a, it returns a value
            normal(a) such that there exists a unit u such that a = u * normal(a). If set to None, an adequate normal
            will be generated if possible, and else a ValueError will be thrown.

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

    if normal is None:
        if domain is ZZ:
            normal = lambda z: z.abs()
        elif is_PolynomialRing(domain) and domain.base().is_field():
            normal = lambda f: f.parent().zero() if f.is_zero() else f.quo_rem(f.lc())[0]
        else:
            raise ValueError("No default implementation for normal found, a value must be provided")

    q = [domain.zero()]
    rho = [__lu(f, normal), __lu(g, normal)]
    r = [normal(f), normal(g)]
    s = [domain.one().quo_rem(rho[0])[0], domain.zero()]
    t = [domain.zero(), domain.one().quo_rem(rho[1])[0]]

    i = 1
    while r[i] != domain.zero():
        q.append(r[i - 1].quo_rem(r[i])[0])
        rho.append(__lu(r[i - 1] - q[i] * r[i], normal))
        r.append((r[i - 1] - q[i] * r[i]).quo_rem(rho[-1])[0])
        s.append((s[i - 1] - q[i] * s[i]).quo_rem(rho[-1])[0])
        t.append((t[i - 1] - q[i] * t[i]).quo_rem(rho[-1])[0])
        i += 1

    return r, s, t, q


def result_to_pretty_string(f, g, rs, ss, ts):
    """
    Creates a nicely formatted string using the result from the algorithm.
    """

    return '\n'.join(map(lambda (r, s, t): '{r} = ({s}) * ({f}) + ({t}) * ({g})'.format(r=r, s=s, t=t, f=f, g=g),
                         zip(rs, ss, ts)))


def main():
    """ Execute the examples. """
    print extended_euclidean_algorithm(ZZ(126), ZZ(35))
    print normalized_extended_euclidean_algorithm(ZZ(126), ZZ(35))

    R = PolynomialRing(QQ, 'x')
    print extended_euclidean_algorithm(R('x^3 - 1'), R('x + 1'))
    print normalized_extended_euclidean_algorithm(R('x^3 - 1'), R('x + 1'))

    S = PolynomialRing(GF(11), 'x')
    print normalized_extended_euclidean_algorithm(S('x'), S(5))


if __name__ == '__main__':
    main()
