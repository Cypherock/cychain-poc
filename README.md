# Decentralised Access Control through Fast TSS
The protocol uses a linearly homomorphic encryption scheme using class groups based on [CL15](https://eprint.iacr.org/2015/047)'s work to compute encrypted signature from the encrypted ECDSA signing key. 

The encryption scheme has multiple advantageous features:

- The message space can be adjusted to any prime, or a power of prime which is suitable for ECDSA circuit computation.
- Because of the correct message space, we don't require additional range proofs.
- The construction of keys for the scheme does require any secret primes.

The scheme is constructed via class groups of unknown order which can be created with a subgroup (F) where the discrete log problem is easy to solve. The other subgroup (H) has an unknown order where the discrete log problem is hard.

## Installation
To compile the code, a C++ compiler and CMake 3.5.1 or later are necessary, and the following libraries are required:
- GMP
- openSSL

On Debian and Ubuntu, the necessary files can be installed with
`apt install g++ libgmp-dev libssl-dev cmake`

### Steps to build:
- Clone the repo and the submodules, `git clone https://github.com/Cypherock/cychain-poc.git --recurse-submodules`
- `cd cychain-poc`
- `mkdir build`
- `cd build`
- `cmake ..`
- `make`

### Execution Instructions
Run the executable `cychain-sig-poc` in `build/` directory as follows:
`./cychain-sig-poc <num-of-users> <num-of-validators> <sec-level>`

Examples:
`./cychain-sig-poc 10 1000 128`
`./cychain-sig-poc 100 1000 128`
`./cychain-sig-poc 100 2000 128`

> \<sec-level\> can only be 128 or 256.

## Cryptography Primitives
### Threshold Additive Homomorphic Encryption
We use the CL-HSMq scheme introduced in [https://eprint.iacr.org/2018/791](https://eprint.iacr.org/2018/791) which has an adjustable message space of an odd prime `q`, but also variant modulo `q^k` and product of primes as analysed in [DJS19](https://dl.acm.org/doi/abs/10.1007/978-3-030-16458-4_20). We define `q` to be the order of `SECP256K1` curve to solve ECDSA circuit in encrypted form. We use the C++ class group implementation [BICYCL](https://eprint.iacr.org/2022/1466) for our PoC which also support CL-HSM2k construction described in [https://eprint.iacr.org/2022/1143](https://eprint.iacr.org/2022/1143) (although not required for us).

The scheme describes the following functions:

#### KeyGen
- Sample private decryption key `sk` randomly in the private key bounds defined by the public parameters.
- Compute public encryption key `pk = h ^ sk`.

#### Encrypt
- Sample a random number `r` in the private key bounds.
- Compute `c1 = h ^ r` and `c2 = f^m * pk^r`.
- Return `(c1, c2)` as the ciphertext.

#### Addition (+)
Given the public key `pk` and encrypted messages `e1` and `e2`, perform the following:

- Parse `(c1, c2) <- e1` and `(c1', c2') <- e2`.
- Compute `c1'' = c1 * c1'` and `c2'' = c2 * c2'`.
- Sample a random number `r` in the private key bounds.
- Return `e = e1 + e2 = (c1'' * h^r, c2'' * pk^r)`.

#### Scalar Multiplication (*)
Given the public key `pk`, encrypted message `e` and a scalar `a`, perform the following:

- Parse `(c1, c2) <- e`.
- Compute `c1' = c1 ^ a` and `c2' = c2 ^ a`.
- Sample a random number `r` in the private key bounds.
- Return `a * e = (c1' * h^r, c2' * pk^r)`.

To support threshold decryption of an encrypted message, we use two more functions, `partDec()` and `aggPartDecs()`.

#### Partial Decryption
Given an encrypted message `e` and a decryption key additive share `sk_i`. (Additive share can be obtained from the shamir integer share by multiplying with the appropriate lagrange coefficient), perform the following:

- Parse `(c1, c2) <- e`.
- Compute the partial decryption `d_i = c1 ^ (-sk_i)` and return this value.

#### Aggregate Partial Decryptions
Given an encrypted message `e` and a set of threshold number of partial decryptions \[`d_i`\], perform the following:

- Parse `(c1, c2) <- e`.
- Compute `M = c2 * product(d_i)`.
- Return `log_f(M)` if `M` is in `F` else return null.

### Zero Knowledge Proofs
We require the following proofs in the protocol:

1. Proof of knowledge of discrete log. (knowledge of `k` given `k.G`)
2. Proof of plaintext knowledge and correct multiplication.
3. Proof of encryption of discrete log. (prove knowledge of `k` and that `k` used in `enc(k)` and `k.G` is same).
4. Proof of correct threshold decryption.

We require these proofs to achieve the identifiable abort property which allows us to detect the malicious or corrupt nodes and prevent them from further performing DoS attacks against the signature generation process. Some proofs might be redundant and the exact proofs need to be finalised.

The paper [https://eprint.iacr.org/2021/205](https://eprint.iacr.org/2021/205) provides the design for 1st and 3rd proofs and the paper [https://eprint.iacr.org/2022/1437](https://eprint.iacr.org/2022/1437) provides the design for the 2nd proof. The design for 4th proof can be constructed from the above proof systems after only a slight modification.

### Integer Secret Sharing
The reconstruction or the usage of shamir shares involves division in the lagrange coefficients which are not supported groups with unknown order. The following protocol by [https://eprint.iacr.org/2022/1437](https://eprint.iacr.org/2022/1437) is used to share secret integers:
- To distribute the secret `s`, create polynomial `F(x) = n! * s + a_{1} * x + ... + a_{t-1} * x^(t - 1)`. (Here, `n!` is multiplied to the integer `s` to prevent the leakage `s mod i` when `F(i)` is shared).
- Share `F(i)` for all `i` in set of participant indices.
- To reconstruct, get `t` or more shares and return `(n!)^2 * s = sum(n! * L_i * F(i))` where `L_i` is the lagrange coefficient. It is multiplied with `n!` to cancel out the denominator.

**Notes**
> The reconstruction returns `(n!)^2 * s` instead of `s`.

> Other coefficients in the polynomial should be chosen from a large enough interval to hide a. The range mentioned by [https://www.ndss-symposium.org/ndss-paper/secure-multiparty-computation-of-threshold-signatures-made-more-efficient/](https://www.ndss-symposium.org/ndss-paper/secure-multiparty-computation-of-threshold-signatures-made-more-efficient/) is `[0, 2^(log2(B) + 1 + 2logt + nlogn + λd)]` where `B` is the secret key bound of class group parameters.

### Chunking
To share class group key shares which are of length `(n^t)*2^(log2(B) + 1 + 2logt + nlogn + λd)` using the encryption scheme with message space `q` where `q` is the order of curve `SECP256K1`, we have to perform chunking of the shares.
- Decompose share `X` in base `q` such that `X = x_k * q^k + x_{k-1} * q^(k-1) + ... + x_1 * q + x_0`. This decomposition can be done by performing repeated mod `q` and integer division.
- Return the chunks `[x_k, x_{k-1}, ..., x_1, x_0]`.

## Distributed Key Generation
We require the distribution of class group decryption key among participants so that they can perform threshold decryption of an encrypted message using their share of the decryption key. 

We refer the 2 round DKG protocol described in [https://www.ndss-symposium.org/ndss-paper/secure-multiparty-computation-of-threshold-signatures-made-more-efficient/](https://www.ndss-symposium.org/ndss-paper/secure-multiparty-computation-of-threshold-signatures-made-more-efficient/) to prevent public key biasing.

> Note: Pederson commitments are used in class groups.

### Round 1
All `n` users perform the following:
- Sample two random secret values `X_i` and `X'_i` in the secretkey bounds of CL cryptosystem. `X'_i`'s shares are used as hiding factors in Pederson commitments.
- Compute `n` shares using the integer secret sharing of the above two secrets. The `j`-th shares are represented as `X_ij` and `X'_ij`.
- For each share `X_ij`, compute it's Pederson commitment using `X'_ij` as `PC_ij`.
- Compute `PK_ij = h^(X_ij * n!)` and encrypt this class group element using `j`-th participant's encryption key to get `C-PK_ij`. Here `h` is the generator of subgroup of unknown order.
- Compute chunks of all shares `X_ij` and encrypt them using `j`-th participant's encryption key.
- Compute the proof `P_ij` to prove the usage of same `X_ij` in the chunks and `PC_ij`.
- Broadcast the proof, encrypted values and Pederson commitments.

At the end, the `n` users peform the following:
- Compute duals for all participant indices.
- Perfrom dual-code verification on all the Pederson commitments `PC_ij` in class groups.
- Verify all the proofs `P_ij`.
- Remove the participants for whose data, the verification fails.

### Round 2
All `n` users perform the following:
- Get all the encrypted chunked shares from `j` parties, and decrypt the chunks.
- Combine the chunks to get `j` shares and add all shares to get one share `x_i` for self.
- Compute `Pub-x_i = h^(x_i * n!)`.
- Combine all `C-PK_ij` to get `C-PK_i` using homomorphic addition.
- Compute proof of knowledge of `x_i` and correct `Pub-x_i` and `C-PK_i`.
- Broadcast the proof, `Pub-x_i` and `C-PK_i`.

At the end, the `n` users perform the following:
- Combine all `C-PK_ij` to get `C-PK_i` using homomorphic addition.
- Verify the proof for `C-PK_i` and `Pub-x_i`.
- Remove the participant `i` if the proof verification fails.
- Return the public encryption key as, `PK = prod(Pub-x_i ^ (n! * L_i))` where `L_i` is lagrange coefficient.

To achieve proactive security, we use the distributed key re-sharing technique desribed in [https://eprint.iacr.org/2021/339](https://eprint.iacr.org/2021/339) but modified for class group keys.

## Definitions
- **User Group** -  Collective of individuals who jointly possess cryptocurrency assets.
- **Policy** - Programmable rules describing conditions that must be fulfilled to generate signatures on each possible message.
- **Policy Contract** - Policy implemented in the form of a smart contract.
- **Validator** - Helps in generating signatures upon successful policy compliance check.
- **Validator Network** - Large network of validators to achieve decentralisation for policy enforcement.

## Validator Network Setup
Validator network is setup by performing the DKG process described above. Finally, a public encryption key of the validator network is generated (`ek_v`) and every validator has a share of decryption key, such that any threshold (pre-defined) number of validators can come online to evaluate the complete decryption key (although never required).

## User Group Setup
User Group first defines the policy and the number of different threshold values (`|T|`) programmed in the policy. The group then performs class group DKG `|T|` number of times with appropriate threshold values (declared in the policy).

After the generation of `|T|` public encryption keys, they are combined with the validator network's public key (`ek_v`) to get `|T|` global encryption keys.

Each user group member then perform the following to generate ECDSA siging key distributively:

- Sample a random number `a_i` less than the order of `SECP256K1` (this will be the additive share of the signing key).
- Compute `a_i.G` and `enc(a_i, ek)` for all `ek` in the global encryption keys.
- Broadcast `a_i.G` and the encrypted values along with the proofs.

Finally, all the encrypted values (for a particular global encryption key) are added together to get the encrypted signing key.
- `enc(X, ek) = sum(enc(a_i, ek))` for all `ek` in global encryption keys (homomorphic addition).

These keys are then broadcasted to the network.

## Access Control Flow
After the setup, the signature over a message can be generated as follows:

- User group creates a signature generation request to the network by providing the message and number of parties available to sign the message.
- The network checks the message and ids of online members against the policy contract for compliance.
- If the compliance check passes, the validator network helps in generating encrypted signature from the encrypted private key and then help in the decryption process.

> Note: There are 3 rounds in the signature generation process, and different set of validators can participate in different rounds, making the scheme more robust and high-throughput.

## Signature Generation
Any `t_p` number of users, if they wish to generate a signature on a message `M` participate in the following three rounds, after successful policy compliance check, using the encrypted signing key `enc(X, ek_tp)` where `ek_tp`'s corresponding decryption key is distributed such that `t_p` number of users + `t_v` number of validators from validator network can decrypt the data.

### Round 1
All `tp` number of users and any set of `t_v` number of validators perform the following:

- Sample a random value `k_i` less than the order of `SECP256K1` curve.
- Compute `K_i = k_i.G` where `G` is the generator point and `.` represents ECC scalar point multiplication.
- Compute `enc(k_i, ek_tp)`.
- Broadcast the values `K_i` and `enc(k_i, ek_tp)`.

At the end, all `K_i` and `enc(k_i, ek_tp)` values are combined to get:

- `K = sum(K_i)` (ECC point addition)
- `r = x-coord(K)`
- `enc(k, ek_tp) = sum(enc(k_i, ek_tp))` (homomorphic addition)

### Round 2
All `tp` number of users and any set of `t_v` number of validators perform the following:

- Sample a random value `p_i` less than the order of `SECP256K1` curve.
- Compute `enc(p_i, ek_tp)`.
- Fetch `enc(k, ek_tp)` and compute `enc(p_i * k, ek_tp) = p_i * enc(k, ek_tp)` (homomorphic scalar multiplication).
- Fetch `enc(X, ek_tp)` and compute `enc(p_i * X, ek_tp) = p_i * enc(X, ek_tp)` (homomorphic scalar multiplication).
- Broadcast the values `enc(p_i, ek_tp)`, `enc(p_i * k, ek_tp)` and `enc(p_i * X, ek_tp)`.

At the end, the following computation is performed:

- `enc(p, ek_tp) = sum(enc(p_i, ek_tp))` (homomorphic addition)
- `enc(p * k, ek_tp) = sum(enc(p_i * k, ek_tp))` (homomorphic addition)
- `enc(p * X, ek_tp) = sum(enc(p_i * X, ek_tp))` (homomorphic addition)
- `enc(z, ek_tp) = (H(M) * enc(p, ek_tp)) + (r * enc(p * X, ek_tp))` (homomorphic addition and scalar multiplication)

### Round 3
All `tp` number of users and any set of `t_v` number of validators perform the following:

- Compute partial decryption `w_i = partDec(enc(p * k, ek_tp))`.
- Compute partial decryption `z_i = partDec(enc(z, ek_tp))`.
- Broadcast the values `w_i` and `z_i`.

At the end, the following computation is performed:

- Aggregate partial decryptions of `enc(p * k, ek_tp)` to get `w = aggPartDecs(w_i)`.
- Aggregate partial decryptions of `enc(z, ek_tp)` to get `z = aggPartDecs(z_i)`.
- Compute `s = z / w` (modulo order of `SECP256K1` curve).

Finally, the `r` and `s` values are combined together to get the signature.

## References
1. [https://eprint.iacr.org/2015/047](https://eprint.iacr.org/2015/047)
2. [https://eprint.iacr.org/2018/791](https://eprint.iacr.org/2018/791)
3. [https://dl.acm.org/doi/abs/10.1007/978-3-030-16458-4_20](https://dl.acm.org/doi/abs/10.1007/978-3-030-16458-4_20)
4. [https://eprint.iacr.org/2022/1466](https://eprint.iacr.org/2022/1466)
5. [https://eprint.iacr.org/2022/1143](https://eprint.iacr.org/2022/1143)
6. [https://eprint.iacr.org/2021/205](https://eprint.iacr.org/2021/205)
7. [https://eprint.iacr.org/2022/1437](https://eprint.iacr.org/2022/1437)
8. [https://www.ndss-symposium.org/ndss-paper/secure-multiparty-computation-of-threshold-signatures-made-more-efficient/](https://www.ndss-symposium.org/ndss-paper/secure-multiparty-computation-of-threshold-signatures-made-more-efficient/)
9. [https://eprint.iacr.org/2021/339](https://eprint.iacr.org/2021/339)
