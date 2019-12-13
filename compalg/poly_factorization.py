"""
Factoring algorithm in Z[x].
"""
from sage.all import *

from extended_euclidean_algorithm import normalized_extended_euclidean_algorithm
from gcd_ufd import gcd_ufd
from finite_field_poly_factorization import berlekamp_poly_factorization
from auxiliary_algorithms import poly_content

# Modern-Computer-Algebra.pdf p.433 seccion 15.1 En Z[x] y Q[x]
# compalg-2017-18.pdf p.55 seccion 2.4.2 Hensel lifting


def __hensel_step(m, f, g, h, s, t):
    """
    Parameters
    ----------
    m: module
    f, g, h, s, t: polynomials in Z[x] such that f = g*h mod m and s*g + t*h = 1 mod m, lc(f) is not a zero divisor
                    modulo m, h is monic, deg(f) = n = deg(g) + deg(h), deg(s) < deg(h) and deg(t) < deg(g)

    Returns
    -------
    Polynomials g', h', s', t' in Z[x] such that f = g'*h' mod m^2 and s'*g' + t'*h' = 1 mod m^2, h' is monic,
    g' = g mod m, h' = h mod m, s' = s mod m, and t' = t mod m, deg(g') = deg(g), deg(h') = deg(h), deg(s') < deg(h'),
    and deg(t') < deg(g').
    """

    R = PolynomialRing(IntegerModRing(m**2), 'x')
    ZZR = PolynomialRing(ZZ, 'x')
    f = R(f)
    g = R(g)
    h = R(h)
    s = R(s)
    t = R(t)

    e = f - g*h
    q, r = (s*e).quo_rem(h)
    g_ = g + t*e + q*g
    h_ = h + r

    b = s*g_ + t*h_ - 1
    c, d = (s*b).quo_rem(h_)
    s_ = s - d
    t_ = t - t*b - c*g_

    return ZZR(g_), ZZR(h_), ZZR(s_), ZZR(t_)


def __multifactor_hensel_lifting(f, p, l, modular_factors):
    """
    Parameters
    ----------
    f: Polynomial in Z[x] of degree n >= 1 such that its leading coefficient is a unit modulo p.
    p: an integer as needed by f and modular_factors.
    l: a natural number as needed by modular_factors.
    modular_factors: monic, nonconstant polynomials f_1, ..., f_r in Z[x] that are pairwise Bezout-coprime modulo p and
                        satisfy f = lc(f) * f_1 * ... * f_r mod p

    Returns
    -------
    Monic polynomials f*_1, ..., f*_r in Z[x] with f = lc(f) * f*_1 * ... * f*_r mod p^l and f*_i = f_i mod p for all i.
    """

    Zp = PolynomialRing(IntegerModRing(p), 'x')
    Zpl = PolynomialRing(IntegerModRing(p**l), 'x')
    ZZR = PolynomialRing(ZZ, 'x')

    r = len(modular_factors)
    if r == 1:
        return [ZZR(Zpl(f / f.leading_coefficient()))]

    k = floor(r / 2)
    d = ceil(log(l, 2))

    modular_factors_p = map(Zp, modular_factors)
    g = ZZR(Zp(f.leading_coefficient()) * prod(modular_factors_p[:k]))
    h = ZZR(prod(modular_factors_p[k:]))

    Fp = PolynomialRing(GF(p), 'x')
    rem, s_list, t_list, _ = normalized_extended_euclidean_algorithm(Fp(g), Fp(h))
    mcd = rem[-2]
    if mcd != ZZR.one():
        raise RuntimeError('Arguments should be coprime, but their mcd is {mcd}'.format(mcd=mcd))
    s = ZZR(s_list[-2])
    t = ZZR(t_list[-2])

    for j in range(d):
        g, h, s, t = __hensel_step(p**(2**j), f, g, h, s, t)

    return __multifactor_hensel_lifting(g, p, l, modular_factors[:k]) \
        + __multifactor_hensel_lifting(h, p, l, modular_factors[k:])


def hensel_lifting_poly_factorization(f):
    """
    Polynomial factorization in Z[x] using Hensel lifting.

    Factors a nonconstant, squarefree, primitive polynomial f in Z[x].

    Parameters
    ----------
    f : A nonconstant, monic polynomial in Z[x]

    Returns
    -------
    The factorization [f_1, ..., f_k] of f.
    """

    domain = f.parent()
    base_domain = domain.base()

    if not base_domain.is_ring() or not base_domain == IntegerRing():
        raise ValueError("The base domain must be the integer ring")

    if f.degree() < 1 or not f.is_squarefree() or not f.is_primitive():
        raise ValueError("f must be a nonconstant, squarefree, primitive polynomial")

    n = f.degree()
    if n == 1:
        return [f]

    A = base_domain(f.norm(Infinity))
    b = f.leading_coefficient()
    B = sqrt(n + 1) * 2**n * A * b
    C = (n + 1)**(2*n) * A**(2*n - 1)
    gamma = ceil(2 * log(C, 2))

    p = 2
    while p <= 2*gamma*log(gamma):
        if b % p != 0:
            Fp = PolynomialRing(GF(p), 'x')
            f_bar = Fp(f)
            if f_bar.is_squarefree():
                break
        p = next_prime(p)

    if p > 2*gamma*log(gamma):  # Should never happen
        raise RuntimeError("Couldn't find such a prime")

    # Modular factorization
    Fp = PolynomialRing(GF(p), 'x')
    f_bar = Fp(f)

    modular_factors = berlekamp_poly_factorization(f_bar, squarefree=True)

    ZZR = PolynomialRing(ZZ, 'x')
    modular_factors = map(ZZR, modular_factors)

    # Hensel lifting
    l = ceil(log(2*B + 1, p))
    modular_factors = __multifactor_hensel_lifting(f, p, l, modular_factors)

    # The set of modular factors still to be treated, the set of factors found, and the polynomial f_ still to be
    # factored.
    Zpl = PolynomialRing(IntegerModRing(p ** l), 'x')
    modular_factors = Set(map(Zpl, modular_factors))
    s = 1
    factors = []
    f_ = f

    # Factor combination
    while 2*s <= len(modular_factors):
        for S in Subsets(modular_factors, s):
            g_ = ZZR(Zpl(b) * prod(S))
            h_ = ZZR(Zpl(b) * prod(modular_factors.difference(S)))

            if g_.norm(1) * h_.norm(1) <= B:
                modular_factors = modular_factors.difference(S)
                factors.append(ZZR(g_ / poly_content(g_)))  # Primitive part
                f_ = ZZR(h_ / poly_content(h_))
                b = f_.leading_coefficient()
                break  # Exit the for loop and continue the while loop

        s += 1

    factors.append(f_)
    return factors


def main():
    """ Execute the examples. """

    R = PolynomialRing(ZZ, 'x')

    f = R('x^4 - 1')
    g = R('x^3 + 2*x^2 - x - 2')
    h = R('x - 2')
    s = R('-2')
    t = R('2*x^2 - 2*x - 1')
    print __hensel_step(5, f, g, h, s, t)  # Expected: (x^3 + 7*x^2 - x - 7, x - 7, 8, -8*x^2 - 12*x - 1)

    # Expected: [x - 1, x - 182, x + 182, x + 1]
    print __multifactor_hensel_lifting(f, 5, 4, [R('x - 1'), R('x - 2'), R('x + 2'), R('x + 1')])

    f2 = R('6*x^4 + 5*x^3 + 15*x^2 + 5*x + 4')
    print hensel_lifting_poly_factorization(f2)  # Expected: [3*x^2 + x + 1, 2*x^2 + x + 4]


if __name__ == '__main__':
    main()
