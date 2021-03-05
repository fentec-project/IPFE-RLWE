# Inner-Product Functional Encryption from Ring-Learning With Errors
Library for software implementation of ring-LWE based inner-production functional encryption (IPFE) 

## Description
Functional encryption (FE) is one of the three components of Computation on Encrypted Data (COED) paradigm. FE allows computing functions on data while maintaining privace of the data. For a given $n$*-dimensional* private integer vector $\bf{X=(x_0, x_1,\cdots, x_{n-1})}$ and a corresponding known weight vector $\bf{Y=(y_0, y_1,\cdots, y_{n-1})}$ and a *secret key* $\bf{sk}$ an IPFE scheme only returns the final inner product of $\bf{x}$ and $\bf{y}$ *i.e* $\bf{F_{sk}(X,Y)=\Sigma_i x_i\cdot y_i}$, whithout revealing anything about vector $\bf{X}$.


Software implementation of inner product functional encryption based on Ring-LWE as a proof of concept. The code supports parameter sets for low, medium and high security, and the implementation is optimized using efficient and constant-time Gaussian sampling, a residue number system to break down modular arithmetic, the Number Theoretic Transform for efficient polynomial multiplication and a pseudo-random number generator based on AES-CTR. The functions for setup, key generation, encryption and decryption can be directly used by any application utilizing inner product functional encryption, taking this implementation as the corresponding cryptographic library.
