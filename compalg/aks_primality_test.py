"""
AKS primality test.
"""
from sage.all import *

from euclidean_algorithm import euclidean_algorithm

import timeit


def __sieve_of_eratosthenes(n):
    """
    Computes a list of prime numbers up to n.

    Complexity: (log n) (log log log n).
    """

    if n < 2:
        return []

    primes = [2]
    candidate = 3

    while candidate < n:
        if all(map(lambda p: candidate % p != 0, primes)):
            primes.append(candidate)

        candidate += 2

    return primes


def __is_perfect_power(n):
    """
    Checks if n can be written as a^b with b > 1.

    Algorithm A from http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.108.458&rep=rep1&type=pdf.
    Complexity (log^3 n) (log log log n).
    """

    max_b = integer_floor(log(n, 2))

    for p in __sieve_of_eratosthenes(max_b):
        x = integer_floor(n ** (1 / float(p)))

        if x ** p == n:
            return True

    return False


def aks_is_prime(n):
    """
    Checks if n is prime.
    """
    # Make sure this is a sage integer
    n = ZZ(n)

    # Input checks
    if n <= 1:
        raise ValueError("The number n must be an integer greater than 2")

    # Step 1: discard perfect powers
    if __is_perfect_power(n):
        return False

    # Step 2: find r
    # o_r(a) = smallest k such as a^k = 1 mod r, when gcd(a, r) = 1
    # We want the smallest r such as o_r(n) > log^2 n
    log_n = log(n, 2)
    log_squared = log_n ** 2

    r = 1
    k = 1
    while k <= log_squared:  # log_n >= 1, therefore k <= log_n and the loop will always run at least once
        r += 1
        power = n % r
        k = 1

        while k <= log_squared:
            if power == 1:
                break     # k < log^2 n but k = o_r(n)

            k += 1
            power = (power * n) % r

        # If the previous loop has ended without finding k = o_r(n), then k > log_squared and the loop will break with
        # the correct value of r. If the previous loop breaks, then we will try another r.

    # Step 3: if 1 < gcd(a, n) < n for any a <= r, output composite
    for a in reversed(IntegerRange(1, r + 1)):
        gcd = euclidean_algorithm(a, n)
        if 1 < gcd < n:
            return False

    # Step 4: if n <= r, output prime
    if n <= r:
        return True

    # Step 5: for a in 1..floor(sqrt(phi(r)) log n), if (x + a)^n != x^n + a mod x^r - 1 with coefficients in Z_n,
    # then output composite.
    # It turns out that calculating phi (Euler's totient function) is expensive, but we can make do with an upper bound
    # with the same time complexity. Actually, we can use r, since phi(r) <= r - sqrt(r) < r.
    Zn = Zmod(n)
    Zn_pol = PolynomialRing(Zn, 'x')
    x = Zn_pol('x')

    for a in list(IntegerRange(integer_floor(sqrt(r) * log_n))):
        first_pol = ((x + Zn_pol(a)) ** n) % (x ** r - Zn_pol(1))
        second_pol = x ** (n % r) + Zn_pol(a)
        if first_pol != second_pol:
            return False

    # Step 6: output prime
    return True


def main():
    """ Execute the examples. """
    prime_list = [2, 11, 101, 1009, 10007, 100003, 1000003, 10000019]
    for p in prime_list:
        start_time = timeit.default_timer()
        if not aks_is_prime(p):
            print 'Something went wrong, {p} is prime'.format(p=p)
        elapsed = timeit.default_timer() - start_time
        print 'Checked primality of {p} in time {t}s'.format(p=p, t=elapsed)

    start, end = 2, 200
    for n in range(start, end):
        print 'Started test for n = {n}'.format(n=n)
        if aks_is_prime(n) != is_prime(n):
            print 'The AKS implementation is wrong for {n}'.format(n=n)

    print 'The implementation is right for primes {s} to {e}'.format(s=start, e=end)


if __name__ == '__main__':
    main()
