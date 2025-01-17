"""
AKS primality test.
"""
from sage.all import *

from a1_euclidean_algorithm import euclidean_algorithm

import timeit
from _auxiliary_algorithms import repeated_square


def miller_rabin_is_prime(n, k=10):
    """
    Checks if n is probably prime, with k tries.

    Returns True if n is probably prime after k tries, False if it is certainly composite.
    """
    # Make sure this is a sage integer
    n = ZZ(n)

    # Input checks
    if n <= 1:
        raise ValueError("The number n must be an integer greater than 2")

    if n == 2 or n == 3:
        return True

    # Write n as 2^r * d + 1, d odd
    d = n - 1
    r = 0
    while d % 2 == 0:
        d = d.quo_rem(2)[0]
        r += 1

    # We try to find a witness k times
    for _ in range(k):
        # Choose a random number, and check if it is a witness
        a = randint(2, n - 2)

        # The test is based on checking a^d mod n, so we calculate it efficiently
        x = repeated_square(a, d, n)

        # For a to be a witness, two conditions must hold:
        # (1) a^d != 1 (mod n)
        # (2) For all 0 <= s <= r - 1, (a^d)^(2^s) != -1 (mod n).

        # We check condition (1) and condition (2) for s = 0
        if x == 1 or x == n - 1:
            continue

        # We check condition 2 for s in [1, n-1]
        is_witness = True
        for _ in range(n - 1):
            x = repeated_square(x, 2, n)

            if x == n - 1:
                is_witness = False
                break   # No point in trying more, it is not a witness

        if is_witness:
            # a is a witness, so the number is composite
            return False

    # If none of the iterations have found a witness, then n is probably prime
    return True


def main():
    """ Execute the examples. """
    miller_rabin_is_prime(5)

    start, end = 2, 200
    for n in range(start, end) + [randint(1000, 100000) for _ in range(100)]:
        probably_prime = miller_rabin_is_prime(n)
        if is_prime(n):
            if probably_prime:
                print 'PROBABLY PRIME for PRIME number {n}'.format(n=n)
            else:
                print '!!!!! COMPOSITE for PRIME number {n}'.format(n=n)
        else:
            if probably_prime:
                print 'PROBABLY PRIME for COMPOSITE number {n}'.format(n=n)
            else:
                print 'COMPOSITE for COMPOSITE number {n}'.format(n=n)


if __name__ == '__main__':
    main()
