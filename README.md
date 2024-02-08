
# What is this repository about
A repository with examples for each tricks used in software polynomial multiplication works.

## Intended role of this repository

## Structure of this repository

- `examples_basic`: Examples of each of the ideas.
    - `BigIntMul.c`: This file demonstrates how to multiply integers via polynomial multiplications.
        - Assumed knowledge: Integer and polynomial multiplications.
        - Reference: Folklore.
    - `DWT.c`: This file demonstrates discrete weighted transform (DWT).
        - Assumed knowledge: Chinese remainder theorem for polynomial rings.
        - Reference:
    - `FNT.c`: This file demonstrates Fermat number transform.
        - Assumed knowledge: Chinese remainder theorem for polynomial rings.
        - Reference:
    - `GT.c`: This file demonstrates Good--Thomas FFT.
        - Assumed knowledge: Multi-variate polynomial rings (minimum); tensor product of associate algebras (recommended).
        - Reference:
    - `Karatsuba.c`: This file demonstrates Karatsuba.
        - Assumed knowledge: Chinese remainder theorem for polynomial rings and evaluation at infinity; module homomorphism (recommended).
        - Reference:
    - `Nussbaumer.c`: This file demonstrates Nussbaumer FFT.
        - Assumed knowledge: Chinese remainder theorem for multi-variate polynomial rings.
        - Reference:
    - `Schoenhage.c`: This file demonstrates Schoenhage FFT.
        - Assumed knowledge: Chinese remainder theorem for multi-variate polynomial rings.
        - Reference:
    - `TMVP.c`: This file demonstrates Toeplitz matrix-vector product built upon Toom-4.
        - Assumed knowledge: Module-theoretic dual of algebra homomorphisms over commutative rings.
        - Reference:

## Software requirements
That's simple. You just need to compile `C` programs.
Since I use `Makefile`, below are the requirements for running the examples:
- `make`
- `gcc`/`clang`

In general, you just need to figure out how to compile `C` programs in your environment.

# References

[Win80] Shmuel Winograd. Arithmetic Complexity of Computations, volume 33. Society for Industrial and Applied Mathematics, 1980. https://epubs.siam.org/doi/10.1137/1.9781611970364.

