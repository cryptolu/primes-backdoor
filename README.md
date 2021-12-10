## About
This repository contains a [SageMath](https://www.sagemath.org) implementation of the prime generation procedure and factorization attack detailed in the paper 

[Factoring Primes to Factor Moduli: Backdooring and Distributed Generation of Semiprimes](https://eprint.iacr.org/2021/1610)

by _Giuseppe Vitto_.


## Functionalities
The script `gen_prime.sage` generates a random prime `p` so that a twist of the Complex Multiplication curve defined over the finite field of size `p` has smooth order with respect to a certain input factor base. Optionally, the `safe_prime` option can be set to `True` to output _safe primes_ instead of just primes.

The script `attack.sage` attempts factorisation of an odd input composite integer `N`, using the attack outlined in the paper. In case `N` has one of its prime factors backdoored as above, this prime will be returned.
