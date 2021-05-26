# Inner-Product Functional Encryption from Ring-Learning With Errors
Library for software implementation of ring-LWE based inner-production functional encryption (IPFE) 

## Description
Functional encryption (FE) is one of the three components of Computation on Encrypted Data (COED) paradigm. FE allows computing functions on data while maintaining privace of the data. For a given $n$*-dimensional* private integer vector $\bf{X=(x_0, x_1,\cdots, x_{n-1})}$ and a corresponding known weight vector $\bf{Y=(y_0, y_1,\cdots, y_{n-1})}$ and a *secret key* $\bf{sk}$ an IPFE scheme only returns the final inner product of $\bf{x}$ and $\bf{y}$ *i.e* $\bf{F_{sk}(X,Y)=\Sigma_i x_i\cdot y_i}$, whithout revealing anything about vector $\bf{X}$.


Software implementation of inner product functional encryption based on Ring-LWE as a proof of concept. The code supports parameter sets for low, medium and high security, and the implementation is optimized using efficient and constant-time Gaussian sampling, a residue number system to break down modular arithmetic, the Number Theoretic Transform for efficient polynomial multiplication and a pseudo-random number generator based on AES-CTR. The functions for setup, key generation, encryption and decryption can be directly used by any application utilizing inner product functional encryption, taking this implementation as the corresponding cryptographic library.

## Requirements

**GMP library** It is necessary for the mutiprecission arithmetic during the decomposition of the polynomials into the residue number system where all operations are performed, as well as to retrieve the final result after decryption. On Linux systems it can be installed as:

> apt-get install libgmp3-dev

**Stack size** The library requires a high amount of stack. In Linux systems the stack size is usually limited to 8192 kb but it can be checked and extended as:

> ulimit -a

> ulimit -s unlimited

## Modules

- **sample** implements the sampling of the public element a

- **gauss** implements the Gaussian sampling used for generating the secrets and error terms

- **aes256ctr** implements a PRF based on AES256 in counter mode

- **randombytes** implements a TRNG based on '/dev/urandom' to sample the seeds that are later used as entropy source in our system

- **ntt** implements the NTT-based polynomial multiplication

- **crt** implements the functions to change domains in the residue number system using the CRT

- **params.h** defines the parameters of the scheme

- **consts.h** defines the constants for the NTT according to the parameters of the scheme

- **arith_rns** implements the modular arithmetic according to the parameters of the system and the chosen moduli

- **rlwe_sife** implements the cryptographic operations using all other modules, and provides the API functions for our IPFE scheme 

## Usage

The sources can be compiled as a library and only the functions defined in 'rlwe_sife.h' need to be accessed to implement all cryptographic operations.

There is an example test program under 'src/test' on how these functions are called.
