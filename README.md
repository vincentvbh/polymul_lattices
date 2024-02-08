
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
        - Reference: [CF94]
    - `FNT.c`: This file demonstrates Fermat number transform.
        - Assumed knowledge: Chinese remainder theorem for polynomial rings.
        - Reference: [AB74]
    - `GT.c`: This file demonstrates Good--Thomas FFT.
        - Assumed knowledge: Multi-variate polynomial rings (minimum); tensor product of associate algebras (recommended).
        - Reference: [Goo58]
    - `Karatsuba.c`: This file demonstrates Karatsuba.
        - Assumed knowledge: Chinese remainder theorem for polynomial rings and evaluation at infinity; module homomorphism (recommended).
        - Reference: [KO62]
    - `Nussbaumer.c`: This file demonstrates Nussbaumer FFT.
        - Assumed knowledge: Chinese remainder theorem for multi-variate polynomial rings.
        - Reference: [Nus80]
    - `Schoenhage.c`: This file demonstrates Schoenhage FFT.
        - Assumed knowledge: Chinese remainder theorem for multi-variate polynomial rings.
        - Reference: [Sch77]
    - `TMVP.c`: This file demonstrates Toeplitz matrix-vector product built upon Toom-4.
        - Assumed knowledge: Module-theoretic dual of algebra homomorphisms over commutative rings.
        - Reference: [Fid73], [Win80]

## Software requirements
That's simple. You just need to compile `C` programs.
Since I use `Makefile`, below are the requirements for running the examples:
- `make`
- `gcc`/`clang`

In general, you just need to figure out how to compile `C` programs in your environment.

# References

[Goo58]
I. J. Good. The Interaction Algorithm and Practical Fourier Analysis. Journal of the Royal Statistical Society: Series B (Methodological), 20(2):361– 372, 1958. https://www.jstor.org/stable/2983896.

[KO62]
A. Karatsuba and Yu. Ofman. Multiplication of many-digital numbers by automatic computers. In Doklady Akademii Nauk, volume 145(2), pages 293–294, 1962. http://cr.yp.to/bib/1963/karatsuba.html.

[Tho63]
Llewellyn Hilleth Thomas. Using a computer to solve problems in physics. Applications of digital computers, pages 44–45, 1963.

[Too63]
Andrei L. Toom. The Complexity of a Scheme of Functional Elements Realizing the Multiplication of Integers. Soviet Mathematics Doklady, 3:714–716, 1963. http://toomandre.com/my-articles/engmat/MULT-E.PDF.

[GS66]
W. M. Gentleman and G. Sande. Fast Fourier Transforms: For Fun and Profit. In Proceedings of the November 7-10, 1966, Fall Joint Computer Conference, AFIPS ’66 (Fall), pages 563–578. Association for Computing
Machinery, 1966. https://doi.org/10.1145/1464291.1464352.

[Rad68]
Charles M. Rader. Discrete Fourier Transforms When the Number of Data Samples Is Prime. Proceedings of the IEEE, 56(6):1107–1108, 1968. https://ieeexplore.ieee.org/document/1448407.

[Fid73]
Charles M. Fiduccia. On the Algebraic Complexity of Matrix Multiplication.
1973. https://cr.yp.to/bib/entries.html#1973/fiduccia-matrix.

[AB74]
Ramesh C. Agarwal and Charles S. Burrus. Fast convolution using Fermat number transforms with applications to digital filtering. IEEE Transactions on Acoustics, Speech, and Signal Processing, 22(2):87–97, 1974. https://ieeexplore.ieee.org/abstract/document/1162555.

[Sch77]
Arnold Schönhage. Schnelle Multiplikation von Polynomen über Körpern der Charakteristik 2. Acta Informatica, 7(4):395–398, 1977. https://link.springer.com/article/10.1007/bf00289470.

[Bru78]
Georg Bruun. z-transform DFT Filters and FFT’s. IEEE Transactions on Acoustics, Speech, and Signal Processing, 26(1):56–63, 1978. https://ieeexplore.ieee.org/document/1163036.

[Win80]
Shmuel Winograd. Arithmetic Complexity of Computations, volume 33. Society for Industrial and Applied Mathematics, 1980. https://epubs.siam.org/doi/10.1137/1.9781611970364.

[Nus80]
Henri J. Nussbaumer. Fast Polynomial Transform Algorithms for Digital Convolution. IEEE Transactions on Acoustics, Speech, and Signal Pro- cessing, 28(2):205–215, 1980. https://ieeexplore.ieee.org/document/1163372.

[CK91]
David G. Cantor and Erich Kaltofen. On Fast Multiplication of Polynomials over Arbitrary Algebras. Acta Informatica, 28(7):693–701, 1991. https://link.springer.com/article/10.1007/BF01178683.

[CF94]
Richard Crandall and Barry Fagin. Discrete Weighted Trans- forms and Large-integer Arithmetic. Mathematics of computa- tion, 62(205):305–324, 1994. https://www.ams.org/journals/mcom/1994-62-205/S0025-5718-1994-1185244-1/?active=current.







