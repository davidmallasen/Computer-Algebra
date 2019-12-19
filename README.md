# Computer-Algebra

Some computer algebra algorithms in SageMath.

## Main

1. Euclid's algorithm for any euclidean domain:
    1. GCD computation.
    2. Traditional extended euclidean algorithm.
    3. Normalized extended euclidean algorithm (for a domain with normal form).
2. Chinese remainder theorem algorithm (compute the inverse).
3. Greatest common divisor in a unique factorization domain.
4. Inverse of an element in a finite field. 
p,f irreducibles in Zp \[x] -> K:=Zp\[x] / (f(x)) |K|=p^(deg f).
5. Irreducibility test of a polynomial in Fq\[x].
6. Discrete logarithm in fields Fq\[x] / (f(x)).
7. Factorization in finite fields:
    1. Factoring algorithm of a polynomial in a finite field parts 1, 2 and 3.
    2. Berlekamp's factorization algorithm in finite fields.
8. Factorization algorithms in Z\[x].
9. AKS primality testing algorithm.

## Extra:

### Groebner basis:

1. Buchberger division algorithm.
2. Ideal membership algorithm.
3. Multivariable division algorithm, with unique remainder

### Primality

1. Miller-Rabin primality test