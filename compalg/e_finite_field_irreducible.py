"""
Irreducibility check in a finite field.
"""
from sage.all import *

from __auxiliary_algorithms import repeated_square, horners_evaluation
from a1_euclidean_algorithm import euclidean_algorithm
from __utils import generate_list, filter_with_mask


def __next_power_of_2(x):
    """
    Returns the smallest r such that x < 2^r.
    """
    return 1 if x == 0 else 1 << (x - 1).bit_length()


def __matrix_blocks(A, m, n):
    """
    Returns a list of four matrix blocks, where the first block is of size (m, n).
    """
    return A[:m, :n], A[:m, n:], A[m:, :n], A[m:, n:]


def __fast_square_matrix_multiplication(A, B):
    """
    Calculates A times B, for A, B square matrices of the same size. It divides both matrices in equal-sized blocks
    (padding with zeroes if needed), then calculates useful matrices of the same size and calls this method recursively
    to multiply them. Finally, the product is recreated from the blocks.
    """

    if A.dimensions() != B.dimensions():
        raise ValueError("The matrices must have the same dimensions")

    if A.dimensions()[0] != A.dimensions()[1]:
        raise ValueError("The matrices are not square")

    n = A.dimensions()[0]

    if n == 1:
        return matrix([[A[0, 0] * B[0, 0]]])

    pad_n = __next_power_of_2(n)

    # Pad the matrices to be of size power of 2
    A = block_diagonal_matrix(A, matrix.zero(pad_n - n), subdivide=False)
    B = block_diagonal_matrix(B, matrix.zero(pad_n - n), subdivide=False)

    m = pad_n / 2
    A11, A12, A21, A22 = __matrix_blocks(A, m, m)
    B11, B12, B21, B22 = __matrix_blocks(B, m, m)

    S1 = A21 + A22
    S2 = S1 - A11
    S3 = A11 - A21
    S4 = A12 - S2

    T1 = B12 - B11
    T2 = B22 - T1
    T3 = B22 - B12
    T4 = T2 - B21

    P1 = __fast_square_matrix_multiplication(A11, B11)
    P2 = __fast_square_matrix_multiplication(A12, B21)
    P3 = __fast_square_matrix_multiplication(S4, B22)
    P4 = __fast_square_matrix_multiplication(A22, T4)
    P5 = __fast_square_matrix_multiplication(S1, T1)
    P6 = __fast_square_matrix_multiplication(S2, T2)
    P7 = __fast_square_matrix_multiplication(S3, T3)

    U1 = P1 + P2
    U2 = P1 + P6
    U3 = U2 + P7
    U4 = U2 + P5
    U5 = U4 + P3
    U6 = U3 - P4
    U7 = U3 + P5

    return block_matrix([[U1, U5], [U6, U7]], subdivide=False)[:n, :n]


def __fast_square_times_arbitrary_matrix_multiplication(A, B):
    """
    A times B, where A is square and B is arbitrary. It uses the fast matrix multiplication algorithm by dividing the
    matrix B into square blocks, with the last block padded by zeroes if necessary.
    """

    if A.dimensions()[0] != A.dimensions()[1]:
        raise ValueError("The first matrix is not square")

    m = A.dimensions()[0]
    n = B.dimensions()[1]
    pad_n = ((n + m - 1) / m) * m     # Integer division, ceil(n/m)

    # Pad B
    B = block_matrix([B, zero_matrix(m, pad_n - n)], nrows=1, subdivide=False)

    blocks_B = [B[:, i:(i+m)] for i in range(0, pad_n, m)]
    blocks_prod = map(lambda block: A * block, blocks_B)
    return block_matrix(blocks_prod, nrows=1, subdivide=False)[:, :n]


def __padded_coefficients(f, pad_to):
    """
    Returns the list of coefficients of a polynomial, adding zeroes if necessary to make the list of size pad_to.
    """
    coefs = f.coefficients(sparse=False)
    return coefs + [f.parent().base().zero()] * (pad_to - len(coefs))


def __fast_modular_composition(f, g, h):
    """
    Computes g(h) mod f. The polynomials must have: deg g, deg h < deg f.

    Parameters
    ----------
    f : a nonzero monic polynomial.
    g : a polynomial.
    h : a polynomial.

    Returns
    -------
    g(h) mod f.
    """

    poly_ring = f.parent()
    n = f.degree()
    m = integer_ceil(sqrt(n))

    g_coefs = g.coefficients(sparse=False)
    gs = [poly_ring(g_coefs[i:(i+m)]) for i in range(0, m * m, m)]

    h_power = poly_ring(1)
    hs = []
    for i in range(m):
        hs.append(h_power)
        h_power *= h
        h_power %= f
    # h_power now has the mth power

    A = matrix(map(lambda h: __padded_coefficients(h, n), hs))
    B = matrix(map(lambda g: __padded_coefficients(g, m), gs))
    BA = __fast_square_times_arbitrary_matrix_multiplication(B, A)

    poly_poly_ring = PolynomialRing(poly_ring, 'y')
    r = poly_poly_ring([poly_ring(list(row)) for row in BA.rows()])
    return horners_evaluation(r, h_power) % f


def __prime_divisors(n):
    """
    Calculates all prime divisors of n naively.
    """
    candidate = 2
    primes = []

    while candidate <= n:
        if n % candidate == 0:
            primes.append(candidate)

            while n % candidate == 0:
                n /= candidate

        candidate += 1

    return primes


def is_irreducible(f):
    """
    Irreducibility check for polynomials in a finite field.

    Parameters
    ----------
    f : the element to be checked. Must be a polynomial with coefficients belonging to a finite field.

    Returns
    -------
    A boolean representing whether f is irreducible or not.
    """
    if f.degree() <= 1:
        return True

    poly_field = f.parent()
    base_field = poly_field.base()

    if not base_field.is_field() or not base_field.is_finite():
        raise ValueError("The base field must be a finite field")

    # The algorithm is based on a corollary stating that a polynomial over a finite field is irreducible iff
    # f | (x^(q^n) - x) and for each prime divisor t of n, gcd(x^(q^(n/t)), f) = 1.

    n = f.degree()
    q = base_field.order()
    x = poly_field(poly_field.variable_name())

    # Step 1: compute x^q mod f, and then a = x^(q^n) mod f. If a != x, then it is reducible
    x_power = repeated_square(x, q, f)

    # The trick here is that we can calculate x^(q^(2^r)) easily for all r, since we can just use fast modular
    # composition to double the exponent r times, starting with x_power.
    x_powers_of_two = generate_list(lambda x: __fast_modular_composition(f, x, x), x_power, n.nbits())

    # Once we have all those easily calculated powers, we can compose the ones that we need: those with a 1 on their
    # binary representation.
    # For example, if n = 5, we take x^q and x^(q^4), we compose them, and then we have (x^q)^(q^4) = x^(q * q^4) =
    # x^(q^(1+4)) = x^(q^5), just what we wanted.
    x_q_n = reduce(lambda a, b: __fast_modular_composition(f, a, b), filter_with_mask(n, x_powers_of_two))

    if x != x_q_n:
        return False    # Reducible

    # Step 2: for all prime divisors...
    for prime in __prime_divisors(n):
        # Step 3: compute b = , using the same trick as in step 1.
        # If gcd(b - x, f) != 1, then it is reducible.
        exp = n / prime
        x_q_exp = reduce(lambda a, b: __fast_modular_composition(f, a, b), filter_with_mask(exp, x_powers_of_two))

        if not euclidean_algorithm(x_q_exp - x, f).is_unit():
            return False    # Reducible

    # Step 4: if we reached this point, it is irreducible
    return True     # Irreducible


def main():
    """ Execute the examples. """
    F1 = GF(7)
    F2 = GF(125)

    R1 = PolynomialRing(F1, 'x')
    R2 = PolynomialRing(F2, 'x')

    # x^4 + 1 reducible over any finite field
    f1 = R1('x^4 + 1')
    f2 = R2('x^4 + 1')

    print is_irreducible(f1)    # False
    print is_irreducible(f2)    # False

    # x^2 - 1 reducible over Z[x] and thus over the finite fields
    f1 = R1('x^2 - 1')
    f2 = R2('x^2 - 1')

    print is_irreducible(f1)    # False
    print is_irreducible(f2)    # False

    for degree in range(1, 10):
        # Some irreducible pols
        f1 = R1.irreducible_element(degree)
        f2 = R2.irreducible_element(degree)

        print is_irreducible(f1)  # True
        print is_irreducible(f2)  # True


if __name__ == '__main__':
    main()
