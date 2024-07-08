#include <string>
#include <iostream>
#include <chrono>
#include <vector>
#include <cmath>
#include <map>
#include <random>
#include <openssl/rand.h>
#include <stdlib.h>
#include <stdint.h>

#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <functional>
#include <atomic>
#include <future>

#include "ThreadPool.h"
#include "CryptoUtils.h"

#include "bicycl.hpp"

extern "C"
{
#include "ecdsa.h"
#include "bignum.h"
#include "curves.h"
#include "bip32.h"
}

using std::string;

typedef struct BIntZKProofGenerationInput {
    std::vector<BICYCL::Mpz> x_l;
    std::vector<BICYCL::Mpz> r_l;
    BICYCL::Mpz r;
    BICYCL::Mpz x_dash;
} BIntZKProofGenerationInput;

typedef struct BIntZKProof {
    BICYCL::Mpz chal;
    BICYCL::QFI C;
    std::vector<BICYCL::QFI> R_l;
    std::vector<BICYCL::QFI> S_l;
    BICYCL::QFI R_0;
    BICYCL::QFI S_0;
    std::vector<BICYCL::Mpz> z1_l;
    std::vector<BICYCL::Mpz> z1_l_reduced;
    std::vector<BICYCL::Mpz> z3_l;
    BICYCL::Mpz z2;
    BICYCL::Mpz z4;
} BIntZKProof;

typedef struct BIntZKProofToProve {
    BICYCL::CL_HSMqk::PublicKey ek;
    BICYCL::QFI PC;
    std::vector<BICYCL::QFI> c_l0;
    std::vector<BICYCL::QFI> c_l1;
    BICYCL::QFI c_0;
    BICYCL::QFI c_1;
} BIntZKProofToProve;

typedef struct GDecCLZKProofGenerationInput {
    BICYCL::CL_HSMqk::SecretKey dk;
    BICYCL::Mpz x;
} GDecCLZKProofGenerationInput;

typedef struct GDecCLZKProof {
    BICYCL::Mpz chal;
    BICYCL::QFI C;
    BICYCL::QFI R;
    BICYCL::QFI S;
    BICYCL::Mpz z1;
    BICYCL::Mpz z2;
} GDecCLZKProof;

typedef struct GDecCLZKProofToProve {
    BICYCL::QFI X;
    BICYCL::QFI c_0;
    BICYCL::QFI c_1;
    BICYCL::CL_HSMqk::PublicKey ek;
} GDecCLZKProofToProve;

class IntegerDualPolynomial {
private:
    unsigned int t; // number of shares required to reconstruct the secret
    unsigned int n;
    BICYCL::Mpz n_fact;
    std::vector<BICYCL::Mpz> coefficients;

    // can make this function parallel
    BICYCL::Mpz factorial(unsigned int n) {
        BICYCL::Mpz result(1UL);
        BICYCL::Mpz temp;

        for (unsigned long i = 2; i <= n; ++i) {
            temp = BICYCL::Mpz (i);
            BICYCL::Mpz::mul(result, result, temp); 
        }

        return result;
    }

    void generateCoefficients(const BICYCL::Mpz &exponent_bound, BICYCL::RandGen &randgen) {
        coefficients.resize(n-t);
        BICYCL::Mpz coeff_bound(exponent_bound);

        for (unsigned int i = 0; i < n-t; ++i) {
            coefficients[i] = randgen.random_mpz(coeff_bound);
        }
    }
public:
    IntegerDualPolynomial(unsigned int t, unsigned int n, const BICYCL::Mpz &exponent_bound, BICYCL::RandGen &randgen)
        : t(t), n(n) {
        
        n_fact = factorial(n);
        generateCoefficients(exponent_bound, randgen);
    }

    IntegerDualPolynomial(unsigned int t, unsigned int n, const std::vector<BICYCL::Mpz> &coefficients)
        : t(t), n(n), coefficients(coefficients) {
        if (coefficients.size() != n-t) {
            throw std::invalid_argument("The number of coefficients must be n-t");
        }

        n_fact = factorial(n);
    }

    const std::vector<BICYCL::Mpz> &get_coefficients() const {
        return coefficients;
    }

    // Method to evaluate the polynomial at a given x value
    BICYCL::Mpz evaluate(const unsigned long x, const std::vector<unsigned int> &P) const {
        BICYCL::Mpz result = coefficients[n-t-1];
        BICYCL::Mpz temp;

        for (int i = n - t - 2; i >= 0; --i) {
            BICYCL::Mpz::mul(temp, result, x); // temp = result * x
            BICYCL::Mpz::add(result, temp, coefficients[i]); // result = temp + coefficients[i]
        }

        BICYCL::Mpz delta_vi = compute_delta_vi(x, P, n_fact);
        BICYCL::Mpz::mul(result, result, delta_vi);

        return result;
    }

    // Static function compute_delta_vi. Ideally, the delta should be computed in the function as n!
    static BICYCL::Mpz compute_delta_vi(int i, const std::vector<unsigned int> &P, const BICYCL::Mpz &delta) {
        // Check if i is in P
        if (std::find(P.begin(), P.end(), i) == P.end()) {
            throw std::invalid_argument("The integer i is not in the list P");
        }

        // Initialize denominator_inverse as Mpz
        BICYCL::Mpz denominator_inverse = delta;
        BICYCL::Mpz temp;

        // Loop through the list of integers
        for (const auto &j : P) {
            if (j != i) {
                // BICYCL::Mpz::mul(numerator, numerator, j);
                BICYCL::Mpz i_minus_j ((signed long) i - j);
                BICYCL::Mpz::divexact(denominator_inverse, denominator_inverse, i_minus_j);
            }
        }

        return denominator_inverse;
    }

    static bool verification(std::vector<BICYCL::Mpz> shares, std::vector<BICYCL::Mpz> duals) {
        if (shares.size() != duals.size()) {
            throw std::invalid_argument("The number of shares and duals must be the same.");
        }

        BICYCL::Mpz result(0UL);
        BICYCL::Mpz temp;

        for (int i = 0; i < shares.size(); ++i) {
            BICYCL::Mpz::mul(temp, shares[i], duals[i]);
            BICYCL::Mpz::add(result, result, temp);
        }

        return result == 0UL;
    }
    
    static bool verification_in_commitment(std::vector<BICYCL::QFI> PCs, 
                                           std::vector<BICYCL::Mpz> duals,
                                           BICYCL::CL_HSMqk &pp) {
        if (PCs.size() != duals.size()) {
            throw std::invalid_argument("The number of shares and duals must be the same.");
        }

        BICYCL::QFI result = pp.Cl_G().one();
        BICYCL::QFI temp;

        for (int i = 0; i < PCs.size(); ++i) {
            pp.Cl_G().nupow(temp, PCs[i], duals[i]);
            pp.Cl_G().nucomp(result, result, temp);
        }

        return result.is_one();
    }
};

class IntegerPolynomial {
private:
    unsigned int t; // number of shares required to reconstruct the secret
    unsigned int n;
    BICYCL::Mpz n_fact;
    std::vector<BICYCL::Mpz> coefficients;

    // can make this function parallel
    static BICYCL::Mpz factorial(unsigned int n) {
        BICYCL::Mpz result(1UL);
        BICYCL::Mpz temp;

        for (unsigned long i = 2; i <= n; ++i) {
            temp = BICYCL::Mpz (i);
            BICYCL::Mpz::mul(result, result, temp); 
        }

        return result;
    }

    void generateCoefficients(const BICYCL::Mpz &exponent_bound, const BICYCL::Mpz &secret, BICYCL::RandGen &randgen) {
        coefficients.resize(t);
        BICYCL::Mpz factorial_n = n_fact;

        BICYCL::Mpz::mul(factorial_n, factorial_n, secret);
        coefficients[0] = factorial_n;

        BICYCL::Mpz coeff_bound(exponent_bound);
        BICYCL::Mpz::mul(coeff_bound, coeff_bound, (unsigned long) 2 * t * t);

        BICYCL::Mpz n_to_n((unsigned long) n);
        BICYCL::Mpz::pow_ui(n_to_n, n_to_n, (unsigned long) n);

        BICYCL::Mpz::mul(coeff_bound, coeff_bound, n_to_n);
        BICYCL::Mpz::mulby2k(coeff_bound, coeff_bound, 42);

        for (unsigned int i = 1; i < t; ++i) {
            coefficients[i] = randgen.random_mpz(coeff_bound);
        }
    }
public:
    // Constructor that takes t, n, exponent bound, and secret
    IntegerPolynomial(unsigned int t, unsigned int n, const BICYCL::Mpz &exponent_bound, const BICYCL::Mpz &secret, BICYCL::RandGen &randgen)
        : t(t), n(n) {
        BICYCL::Mpz factorial_n = factorial(n);
        this->n_fact = factorial_n;
        generateCoefficients(exponent_bound, secret, randgen);
    }

    // Constructor that takes only t, n and exponent bound
    IntegerPolynomial(unsigned int t, unsigned int n, const BICYCL::Mpz &exponent_bound, BICYCL::RandGen &randgen)
        : t(t), n(n) {
        BICYCL::Mpz secret = randgen.random_mpz(exponent_bound);
        BICYCL::Mpz factorial_n = factorial(n);
        this->n_fact = factorial_n;
        generateCoefficients(exponent_bound, secret, randgen);
    }
        
    // Constructor that takes t, n, and a vector of coefficients
    IntegerPolynomial(unsigned int t, unsigned int n, const std::vector<BICYCL::Mpz> &coefficients)
        : t(t), n(n), coefficients(coefficients) {
        if (coefficients.size() != t) {
            throw std::invalid_argument("The number of coefficients must be t");
        }
        // multiply n_fact with the secret
        BICYCL::Mpz factorial_n = factorial(n);
        this->n_fact = factorial_n;
        BICYCL::Mpz::mul(this->coefficients[0], this->coefficients[0], factorial_n);
    }

    // TODO - remove this get coefficients
    const std::vector<BICYCL::Mpz> &get_coefficients() const {
        return coefficients;
    }
        
    // Method to evaluate the polynomial at a given x value
    BICYCL::Mpz evaluate(const unsigned long x) const {
        BICYCL::Mpz result = coefficients[t-1];
        BICYCL::Mpz temp;

        for (int i = t - 2; i >= 0; --i) {
            BICYCL::Mpz::mul(temp, result, x); // temp = result * x
            BICYCL::Mpz::add(result, temp, coefficients[i]); // result = temp + coefficients[i]
        }

        return result;
    }

    // Static function compute_delta_lambda. Ideally, the delta should be computed in the function as n!
    static BICYCL::Mpz compute_delta_lambda(int i, const std::vector<unsigned int> &P, const BICYCL::Mpz &delta) {
        // Check if i is in P
        if (std::find(P.begin(), P.end(), i) == P.end()) {
            throw std::invalid_argument("The integer i is not in the list P");
        }

        // Initialize numerator and denominator_inverse as Mpz
        BICYCL::Mpz numerator(1UL);
        BICYCL::Mpz denominator_inverse = delta;
        BICYCL::Mpz temp;

        // Loop through the list of integers
        for (const auto &j : P) {
            if (j != i) {
                BICYCL::Mpz::mul(numerator, numerator, j);
                BICYCL::Mpz j_minus_i ((signed long) j - i);
                BICYCL::Mpz::divexact(denominator_inverse, denominator_inverse, j_minus_i);
            }
        }

        // Multiply both numerator and denominator_inverse
        BICYCL::Mpz::mul(numerator, numerator, denominator_inverse);

        return numerator;
    }

    // static function that reconstructs secret from shares
    static BICYCL::Mpz reconstruct(unsigned int t, unsigned int n, 
                                   std::vector<BICYCL::Mpz> shares, 
                                   std::vector<unsigned int> indices) {
        if (shares.size() < t || shares.size() > n || shares.size() != indices.size()) {
            throw std::invalid_argument("The number of shares must be >= t and <= n.");
        }

        BICYCL::Mpz result(0UL);
        BICYCL::Mpz delta = factorial(n);

        for (int i = 0; i < shares.size(); ++i) {
            BICYCL::Mpz product(1UL);
            BICYCL::Mpz lambda = compute_delta_lambda(indices[i], indices, delta);

            BICYCL::Mpz::mul(product, product, shares[i]);
            BICYCL::Mpz::mul(product, product, lambda);

            BICYCL::Mpz::add(result, result, product);
        }

        BICYCL::Mpz::divexact(result, result, delta);
        BICYCL::Mpz::divexact(result, result, delta);

        return result;
    }
};

class ClassGroupCommitmentBuilder {
public:
    ClassGroupCommitmentBuilder(const BICYCL::QFI &h, const BICYCL::CL_HSMqk &pp) : h_(h), pp_(pp) {
        precompute();
    }

    void compute_commitment(BICYCL::QFI &r, const BICYCL::Mpz &secret, BICYCL::Mpz &random) {
        BICYCL::QFI base_power;
        pp_.Cl_G().nupow (base_power, h_, secret, d_, e_, h_e_precomp_, h_d_precomp_, h_de_precomp_);

        BICYCL::QFI random_power;
        pp_.power_of_h(random_power, random);

        pp_.Cl_G().nucomp(r, base_power, random_power);
    }

    const BICYCL::QFI &h() const {
        return h_;
    }

private:
    BICYCL::QFI h_;
    BICYCL::CL_HSMqk pp_;
    size_t d_;
    size_t e_;
    BICYCL::QFI h_de_precomp_;
    BICYCL::QFI h_d_precomp_;
    BICYCL::QFI h_e_precomp_;

    void precompute() {
        d_ = (pp_.encrypt_randomness_bound().nbits() + 1) / 2;
        e_ = d_ / 2 + 1;
        h_de_precomp_ = h_;

        for (size_t i = 0; i < d_ + e_; i++) {
            if (i == e_)
                h_e_precomp_ = h_de_precomp_;
            if (i == d_)
                h_d_precomp_ = h_de_precomp_;
            pp_.Cl_G().nudupl(h_de_precomp_, h_de_precomp_);
        }
    }
};

// can make this function parallel
BICYCL::Mpz factorial(unsigned int n) {
    BICYCL::Mpz result(1UL);
    BICYCL::Mpz temp;

    for (unsigned long i = 2; i <= n; ++i) {
        temp = BICYCL::Mpz (i);
        BICYCL::Mpz::mul(result, result, temp); 
    }

    return result;
}

// Function to divide the big number into smaller chunks less than q
std::vector<BICYCL::Mpz> to_q_ary(const BICYCL::Mpz &number, const BICYCL::Mpz &q) {
    std::vector<BICYCL::Mpz> result;
    BICYCL::Mpz current = number;
    BICYCL::Mpz chunk;

    int i = 0;

    while (current > 0UL) {
        BICYCL::Mpz::mod(chunk, current, q); // chunk = current % q
        result.push_back(chunk);

        BICYCL::Mpz::sub(current, current, chunk); // current = current - chunk
        BICYCL::Mpz::divexact(current, current, q); // current = current / q
    }

    return result;
}

// Function to reconstruct the big number from chunks
BICYCL::Mpz from_q_ary(const std::vector<BICYCL::Mpz> &chunks, const BICYCL::Mpz &q) {
    BICYCL::Mpz result(0UL);
    BICYCL::Mpz base(1UL);
    BICYCL::Mpz temp;

    for (const auto &chunk : chunks) {
        BICYCL::Mpz::mul(temp, chunk, base); // temp = chunk * base
        BICYCL::Mpz::add(result, result, temp); // result = result + temp
        BICYCL::Mpz::mul(base, base, q); // base = base * q
    }

    return result;
}

void group_element_encryption(BICYCL::QFI &c1, BICYCL::QFI &c2, 
                              const BICYCL::CL_HSMqk::PublicKey &ek, 
                              const BICYCL::QFI &message, 
                              const BICYCL::Mpz &random, 
                              const BICYCL::CL_HSMqk &pp) {
    pp.power_of_h(c1, random);
    ek.exponentiation(pp, c2, random);
    pp.Cl_Delta().nucomp (c2, c2, message);
}

void group_element_decryption(BICYCL::QFI &message, 
                              const BICYCL::QFI &c1, 
                              const BICYCL::QFI &c2, 
                              const BICYCL::CL_HSMqk::SecretKey &sk, 
                              const BICYCL::CL_HSMqk &pp) {
    pp.Cl_G().nupow (message, c1, sk);
    pp.Cl_Delta().nucompinv (message, c2, message);
}

BICYCL::Mpz generate_BInt_fiat_shamir_hash() {
    BICYCL::Mpz hash;
    return hash;
}

BIntZKProof generate_BInt_proof(BIntZKProofGenerationInput input, 
                                BIntZKProofToProve to_prove, 
                                BICYCL::Mpz &delta,
                                BICYCL::CL_HSMqk &pp, 
                                ClassGroupCommitmentBuilder &cgcb,
                                BICYCL::RandGen &randgen) {
    BIntZKProof proof;

    std::vector<BICYCL::Mpz> a1_l;
    std::vector<BICYCL::Mpz> a3_l;

    BICYCL::Mpz rand_int_bound = pp.encrypt_randomness_bound();
    BICYCL::Mpz::mulby2k(rand_int_bound, rand_int_bound, 128 + 42);

    BICYCL::Mpz a2 = randgen.random_mpz(rand_int_bound);
    BICYCL::Mpz a4 = randgen.random_mpz(rand_int_bound);

    for (size_t i = 0; i < to_prove.c_l0.size(); ++i) {
        a1_l.push_back(randgen.random_mpz(pp.q()));
        a3_l.push_back(randgen.random_mpz(rand_int_bound));
    }

    BICYCL::Mpz a1 = from_q_ary(a1_l, pp.q());
    cgcb.compute_commitment(proof.C, a1, a2);

    for (size_t i = 0; i < to_prove.c_l0.size(); ++i) {
        BICYCL::CL_HSMqk::CipherText ct = pp.encrypt(to_prove.ek, 
                                                     BICYCL::CL_HSMqk::ClearText(pp, a1_l[i]), 
                                                     a3_l[i]);
        proof.R_l.push_back(ct.c1());
        proof.S_l.push_back(ct.c2());
    }

    BICYCL::Mpz delta_a1 = a1;
    BICYCL::Mpz::mul(delta_a1, delta_a1, delta);
    BICYCL::QFI h_to_pow_of_delta_a1;
    pp.power_of_h(h_to_pow_of_delta_a1, delta_a1);

    group_element_encryption(proof.R_0, proof.S_0, to_prove.ek, h_to_pow_of_delta_a1, a4, pp);

    BICYCL::Mpz chal(123UL);
    proof.chal = chal;

    for (size_t i = 0; i < to_prove.c_l0.size(); ++i) {
        BICYCL::Mpz z1_l, z3_l;

        BICYCL::Mpz::mul(z1_l, chal, input.x_l[i]);
        BICYCL::Mpz::add(z1_l, z1_l, a1_l[i]);
        proof.z1_l.push_back(z1_l);

        BICYCL::Mpz::mod(z1_l, z1_l, pp.q());

        BICYCL::Mpz::mul(z3_l, chal, input.r_l[i]);
        BICYCL::Mpz::add(z3_l, z3_l, a3_l[i]);

        proof.z1_l_reduced.push_back(z1_l);
        proof.z3_l.push_back(z3_l);
    }

    BICYCL::Mpz::mul(proof.z2, chal, input.x_dash);
    BICYCL::Mpz::add(proof.z2, proof.z2, a2);

    BICYCL::Mpz::mul(proof.z4, chal, input.r);
    BICYCL::Mpz::add(proof.z4, proof.z4, a4);

    return proof;
}

bool verify_BInt_proof(BIntZKProof proof,
                       BIntZKProofToProve to_prove,
                       BICYCL::Mpz delta,
                       BICYCL::CL_HSMqk &pp,
                       ClassGroupCommitmentBuilder &cgcb) {
    BICYCL::Mpz chal(123UL);

    if (proof.chal != chal) {
        return false;
    }

    BICYCL::Mpz check_bound = pp.encrypt_randomness_bound();
    BICYCL::Mpz::mulby2k(check_bound, check_bound, 128);

    BICYCL::Mpz temp(1UL);
    BICYCL::Mpz::mulby2k(temp, temp, 42);
    BICYCL::Mpz::add(temp, temp, 1UL);

    BICYCL::Mpz::mul(check_bound, check_bound, temp);

    if (proof.z2 < 0UL || proof.z2 > check_bound) {
        return false;
    }

    if (proof.z4 < 0UL || proof.z4 > check_bound) {
        return false;
    }

    for (size_t i = 0; i < proof.z1_l.size(); ++i) {
        if (proof.z1_l_reduced.at(i) < 0UL || proof.z1_l_reduced.at(i) >= pp.q()) {
            return false;
        }

        if (proof.z3_l.at(i) < 0UL || proof.z3_l.at(i) > check_bound) {
            return false;
        }

        BICYCL::Mpz z1_l_reduced;
        BICYCL::Mpz::mod(z1_l_reduced, proof.z1_l.at(i), pp.q());

        if (z1_l_reduced != proof.z1_l_reduced.at(i)) {
            return false;
        }

        BICYCL::CL_HSMqk::CipherText ct = pp.encrypt(to_prove.ek, 
                                                     BICYCL::CL_HSMqk::ClearText(pp, proof.z1_l_reduced[i]), 
                                                     proof.z3_l[i]);
        BICYCL::QFI Rl_times_c_l0_chal = to_prove.c_l0[i];
        BICYCL::QFI Sl_times_c_l1_chal = to_prove.c_l1[i];

        pp.Cl_G().nupow(Rl_times_c_l0_chal, Rl_times_c_l0_chal, chal);
        pp.Cl_Delta().nupow(Sl_times_c_l1_chal, Sl_times_c_l1_chal, chal);

        pp.Cl_G().nucomp(Rl_times_c_l0_chal, Rl_times_c_l0_chal, proof.R_l[i]);
        pp.Cl_Delta().nucomp(Sl_times_c_l1_chal, Sl_times_c_l1_chal, proof.S_l[i]);

        if (!(Rl_times_c_l0_chal == ct.c1()) || !(Sl_times_c_l1_chal == ct.c2())) {
            return false;
        }
    }

    BICYCL::Mpz z1 = from_q_ary(proof.z1_l, pp.q()); 
    BICYCL::Mpz delta_z1 = z1;
    BICYCL::Mpz::mul(delta_z1, delta_z1, delta);
    BICYCL::QFI h_to_pow_of_delta_z1;
    pp.power_of_h(h_to_pow_of_delta_z1, delta_z1);

    BICYCL::QFI R_0_times_c_0_chal = to_prove.c_0;
    BICYCL::QFI S_0_times_c_1_chal = to_prove.c_1;

    BICYCL::QFI lhs_c0;
    BICYCL::QFI lhs_c1;

    group_element_encryption(lhs_c0, lhs_c1, to_prove.ek, h_to_pow_of_delta_z1, proof.z4, pp);

    pp.Cl_G().nupow(R_0_times_c_0_chal, R_0_times_c_0_chal, chal);
    pp.Cl_Delta().nupow(S_0_times_c_1_chal, S_0_times_c_1_chal, chal);

    pp.Cl_G().nucomp(R_0_times_c_0_chal, R_0_times_c_0_chal, proof.R_0);
    pp.Cl_Delta().nucomp(S_0_times_c_1_chal, S_0_times_c_1_chal, proof.S_0);

    if (!(R_0_times_c_0_chal == lhs_c0) || !(S_0_times_c_1_chal == lhs_c1)) {
        return false;
    }

    BICYCL::QFI lhs_pc;
    cgcb.compute_commitment(lhs_pc, z1, proof.z2);

    BICYCL::QFI C_times_PC_chal = to_prove.PC;
    pp.Cl_G().nupow(C_times_PC_chal, C_times_PC_chal, chal);
    pp.Cl_G().nucomp(C_times_PC_chal, C_times_PC_chal, proof.C);

    if (!(C_times_PC_chal == lhs_pc)) {
        return false;
    }

    return true;
}

BICYCL::Mpz generate_GDecCL_fiat_shamir_hash() {
    BICYCL::Mpz hash;
    return hash;
}


GDecCLZKProof generate_GDecCL_proof(GDecCLZKProofGenerationInput input,
                                    GDecCLZKProofToProve to_prove,
                                    BICYCL::Mpz &delta,
                                    BICYCL::CL_HSMqk &pp,
                                    BICYCL::RandGen &randgen) {
    GDecCLZKProof proof;

    BICYCL::Mpz rand_int_bound = pp.encrypt_randomness_bound();
    BICYCL::Mpz::mulby2k(rand_int_bound, rand_int_bound, 128 + 42);

    BICYCL::Mpz a1 = randgen.random_mpz(rand_int_bound);
    BICYCL::Mpz a2 = randgen.random_mpz(rand_int_bound);

    BICYCL::Mpz a2_times_delta = a2;
    BICYCL::Mpz::mul(a2_times_delta, a2_times_delta, delta);

    pp.power_of_h(proof.C, a2_times_delta);
    pp.power_of_h(proof.R, a1);

    proof.S = to_prove.c_0;
    pp.Cl_G().nupow(proof.S, proof.S, a1);
    pp.Cl_Delta().nucomp(proof.S, proof.S, proof.C);

    BICYCL::Mpz chal(123UL);
    proof.chal = chal;

    proof.z1 = input.dk;
    BICYCL::Mpz::mul(proof.z1, proof.z1, chal);
    BICYCL::Mpz::add(proof.z1, proof.z1, a1);

    proof.z2 = input.x;
    BICYCL::Mpz::mul(proof.z2, proof.z2, chal);
    BICYCL::Mpz::add(proof.z2, proof.z2, a2);

    return proof;
}

bool verify_GDecCL_proof(GDecCLZKProof proof,
                         GDecCLZKProofToProve to_prove,
                         BICYCL::Mpz &delta,
                         BICYCL::CL_HSMqk &pp) {
    BICYCL::Mpz chal(123UL);

    if (proof.chal != chal) {
        return false;
    }

    BICYCL::Mpz check_bound = pp.encrypt_randomness_bound();
    BICYCL::Mpz::mulby2k(check_bound, check_bound, 128);

    BICYCL::Mpz temp(1UL);
    BICYCL::Mpz::mulby2k(temp, temp, 42);
    BICYCL::Mpz::add(temp, temp, 1UL);

    BICYCL::Mpz::mul(check_bound, check_bound, temp);

    if (proof.z1 < 0UL || proof.z1 > check_bound) {
        return false;
    }

    if (proof.z2 < 0UL || proof.z2 > check_bound) {
        return false;
    }

    BICYCL::QFI C_times_X_chal = to_prove.X;
    pp.Cl_G().nupow(C_times_X_chal, C_times_X_chal, chal);
    pp.Cl_G().nucomp(C_times_X_chal, C_times_X_chal, proof.C);

    BICYCL::QFI S_times_c1_chal = to_prove.c_1;
    pp.Cl_Delta().nupow(S_times_c1_chal, S_times_c1_chal, chal);
    pp.Cl_Delta().nucomp(S_times_c1_chal, S_times_c1_chal, proof.S);

    BICYCL::QFI R_times_ek_chal;
    to_prove.ek.exponentiation(pp, R_times_ek_chal, chal);
    pp.Cl_G().nucomp(R_times_ek_chal, R_times_ek_chal, proof.R);

    BICYCL::Mpz delta_z2 = proof.z2;
    BICYCL::Mpz::mul(delta_z2, delta_z2, delta);

    BICYCL::QFI h_to_pow_of_delta_z2;
    pp.power_of_h(h_to_pow_of_delta_z2, delta_z2);

    if (!(h_to_pow_of_delta_z2 == C_times_X_chal)) {
        return false;
    }

    BICYCL::QFI h_delta_z2_times_c0_z1 = to_prove.c_0;
    pp.Cl_G().nupow(h_delta_z2_times_c0_z1, h_delta_z2_times_c0_z1, proof.z1);
    pp.Cl_G().nucomp(h_delta_z2_times_c0_z1, h_delta_z2_times_c0_z1, h_to_pow_of_delta_z2);

    if (!(h_delta_z2_times_c0_z1 == S_times_c1_chal)) {
        return false;
    }

    BICYCL::QFI h_to_pow_of_z1;
    pp.power_of_h(h_to_pow_of_z1, proof.z1);

    if (!(h_to_pow_of_z1 == R_times_ek_chal)) {
        return false;
    }

    return true;
}

void unit_test_q_ary(BICYCL::CL_HSMqk &pp) {
    BICYCL::Mpz number("1");

    // Divide number into chunks
    std::vector<BICYCL::Mpz> chunks = to_q_ary(number, pp.q());
    
    // Reconstruct the number from chunks
    BICYCL::Mpz reconstructed_number = from_q_ary(chunks, pp.q());

    if (number == reconstructed_number) {
        std::cout << "Unit test for q_ary passed.\n";
    } else {
        std::cout << "Unit test for q_ary failed.\n";
    }
}

void unit_test_shares_and_duals(BICYCL::CL_HSMqk &pp, BICYCL::RandGen &randgen) {
    unsigned int t = 20;
    unsigned int n = 30;

    BICYCL::Mpz secret(1232346886UL);
    IntegerPolynomial ip = IntegerPolynomial(t, n, pp.secretkey_bound(), secret, randgen);

    std::vector<BICYCL::Mpz> shares;
    std::vector<unsigned int> indices;

    for (int i = 1; i <= n; ++i) {
        shares.push_back(ip.evaluate(i));
        indices.push_back(i);
    }

    BICYCL::Mpz reconstructed = IntegerPolynomial::reconstruct(t, n, shares, indices);

    if (secret == reconstructed) {
        std::cout << "Unit test for share reconstruction passed.\n";
    } else {
        std::cout << "Unit test for share reconstruction failed.\n";
    }

    IntegerDualPolynomial idp = IntegerDualPolynomial(t, n, pp.secretkey_bound(), randgen);

    std::vector<BICYCL::Mpz> duals;

    for (int i = 1; i <= n; ++i) {
        duals.push_back(idp.evaluate(i, indices));
    }

    if (IntegerDualPolynomial::verification(shares, duals) == true) {
        std::cout << "Unit test for dual verification passed.\n";
    } else {
        std::cout << "Unit test for dual verification failed.\n";
    }
}

void unit_test_group_element_encryption(BICYCL::CL_HSMqk &pp, BICYCL::RandGen &randgen) {
    BICYCL::QFI h3; 
    pp.power_of_h(h3, randgen.random_mpz(pp.encrypt_randomness_bound()));

    BICYCL::CL_HSMqk::SecretKey sk = pp.keygen(randgen);
    BICYCL::CL_HSMqk::PublicKey pk = pp.keygen(sk);

    BICYCL::QFI c1, c2;
    group_element_encryption(c1, c2, pk, h3, randgen.random_mpz(pp.encrypt_randomness_bound()), pp);

    BICYCL::QFI message;
    group_element_decryption(message, c1, c2, sk, pp);

    // check if decrypted message is equal to h3
    if (message == h3) {
        std::cout << "Unit test for group element encryption passed.\n";
    } else {
        std::cout << "Unit test for group element encryption failed.\n";
    }
}

void unit_test_class_group_committment(BICYCL::CL_HSMqk pp, BICYCL::RandGen randgen) {
    BICYCL::QFI h3; 
    pp.power_of_h(h3, randgen.random_mpz(pp.encrypt_randomness_bound()));

    BICYCL::Mpz random = randgen.random_mpz(pp.encrypt_randomness_bound());

    BICYCL::QFI commitment1;
    BICYCL::QFI commitment2;

    BICYCL::Mpz secret(123456UL);
    ClassGroupCommitmentBuilder cgc(h3, pp);
    cgc.compute_commitment(commitment1, secret, random);
    cgc.compute_commitment(commitment2, secret, random);

    // check if they are equal
    if (commitment1 == commitment2) {
        std::cout << "Unit test for class group commitment passed.\n";
    } else {
        std::cout << "Unit test for class group commitment failed.\n";
    }
}

void unit_test_BIntZK(BICYCL::CL_HSMqk &pp, BICYCL::RandGen &randgen) {
    unsigned int n = 1000;

    BICYCL::Mpz X = randgen.random_mpz(pp.encrypt_randomness_bound());
    BICYCL::Mpz X_dash = randgen.random_mpz(pp.encrypt_randomness_bound());

    BICYCL::QFI h;
    pp.power_of_h(h, randgen.random_mpz(pp.encrypt_randomness_bound()));

    ClassGroupCommitmentBuilder cgcb(h, pp);
    BICYCL::QFI PC;

    cgcb.compute_commitment(PC, X, X_dash);

    BICYCL::Mpz delta = factorial(n);
    BICYCL::Mpz delta_times_X = delta;
    BICYCL::Mpz::mul(delta_times_X, delta_times_X, X);

    BICYCL::QFI h_to_pow_of_delta_X;
    pp.power_of_h(h_to_pow_of_delta_X, delta_times_X);

    BICYCL::CL_HSMqk::SecretKey sk = pp.keygen(randgen);
    BICYCL::CL_HSMqk::PublicKey ek = pp.keygen(sk);

    BICYCL::QFI c0, c1;
    BICYCL::Mpz r = randgen.random_mpz(pp.encrypt_randomness_bound());
    group_element_encryption(c0, c1, ek, h_to_pow_of_delta_X, r, pp);

    std::vector<BICYCL::Mpz> x_l = to_q_ary(X, pp.q());
    std::vector<BICYCL::Mpz> r_l;

    std::vector<BICYCL::QFI> c_l0;
    std::vector<BICYCL::QFI> c_l1;

    for (size_t i = 0; i < x_l.size(); ++i) {
        r_l.push_back(randgen.random_mpz(pp.encrypt_randomness_bound()));
        BICYCL::CL_HSMqk::CipherText ct = pp.encrypt(ek, BICYCL::CL_HSMqk::ClearText(pp, x_l[i]), r_l[i]);

        c_l0.push_back(ct.c1());
        c_l1.push_back(ct.c2());
    }
        
    BIntZKProofGenerationInput input {x_l, r_l, r, X_dash};
    BIntZKProofToProve to_prove {ek, PC, c_l0, c_l1, c0, c1};

    BIntZKProof proof = generate_BInt_proof(input, to_prove, delta, pp, cgcb, randgen);

    bool result = verify_BInt_proof(proof, to_prove, delta, pp, cgcb);

    if (result == true) {
        std::cout << "Unit test for BInt ZKProof verification passed.\n";
    } else {
        std::cout << "Unit test for BInt ZKProof verification failed.\n";
    }
}

void unit_test_GDecCLZK(BICYCL::CL_HSMqk &pp, BICYCL::RandGen &randgen) {
    unsigned int n = 1000;
    BICYCL::Mpz x = randgen.random_mpz(pp.encrypt_randomness_bound());

    BICYCL::Mpz delta = factorial(n);
    BICYCL::Mpz delta_times_x = delta;
    BICYCL::Mpz::mul(delta_times_x, delta_times_x, x);

    BICYCL::QFI h_to_pow_of_delta_x;
    pp.power_of_h(h_to_pow_of_delta_x, delta_times_x);

    BICYCL::CL_HSMqk::SecretKey sk = pp.keygen(randgen);
    BICYCL::CL_HSMqk::PublicKey ek = pp.keygen(sk);

    BICYCL::QFI c0, c1;
    BICYCL::Mpz r = randgen.random_mpz(pp.encrypt_randomness_bound());
    group_element_encryption(c0, c1, ek, h_to_pow_of_delta_x, r, pp);

    GDecCLZKProofGenerationInput input {sk, x};
    GDecCLZKProofToProve to_prove {h_to_pow_of_delta_x, c0, c1, ek};

    GDecCLZKProof proof = generate_GDecCL_proof(input, to_prove, delta, pp, randgen);

    bool result = verify_GDecCL_proof(proof, to_prove, delta, pp);

    if (result == true) {
        std::cout << "Unit test for GDec-Cl ZKProof verification passed.\n";
    } else {
        std::cout << "Unit test for GDec-Cl ZKProof verification failed.\n";
    }
}

int main() {
    size_t k = 1;
    const ecdsa_curve *curve = get_curve_by_name(SECP256K1_NAME)->params;

    // set seed as time
    BICYCL::Mpz seed = BICYCL::Mpz((unsigned long) time(NULL));
    BICYCL::RandGen randgen(seed);

    // Setting up the the public parameters where message space is the curve order of SECP256K1
    BICYCL::Mpz q(0UL);
    q = string("115792089237316195423570985008687907852837564279074904382605163141518161494337");

    BICYCL::Mpz p(0UL);

    p = string("43193454325736827482201750544140481404008634846532639450547372755026427405991826700723249528255049592325421795981143653914129884011306525120950744025035685855512695631736186909872842154602977956593529422235466196411674948283663780511726036070197880528511138509668014828860898402959171280613034578545771859496996085810790266951981368946435378708249190772155224039947314675049652489452429836293254671938998647344790956467563810567741334084883797963767053599118688017207673143");

    std::cout << "Setting up class group...\n";
    BICYCL::CL_HSMqk pp = BICYCL::CL_HSMqk(q, k, p);
    std::cout << "Class group successfully setup.\n";

    unit_test_q_ary(pp);
    unit_test_shares_and_duals(pp, randgen);
    unit_test_group_element_encryption(pp, randgen);
    unit_test_class_group_committment(pp, randgen);
    unit_test_BIntZK(pp, randgen);
    unit_test_GDecCLZK(pp, randgen);

    std::cout << "\nPerforming DKG-CL mock\n";

    // NOTE: Users are indexed from 1 to n (inclusive).
    unsigned int t = 3;
    unsigned int n = 5;

    // delta = n!
    BICYCL::Mpz delta = factorial(n);

    std::vector<BICYCL::CL_HSMqk::SecretKey> users_secret_keys;
    std::vector<BICYCL::CL_HSMqk::PublicKey> users_public_keys;

    std::vector<unsigned int> indices;

    // users create their own key pairs and broadcast the public keys.
    std::cout << "Users individually create their class group key pairs and broadcast the public keys.\n";
    for (size_t i = 0; i < n; ++i) {
        users_secret_keys.push_back(pp.keygen(randgen));
        users_public_keys.push_back(pp.keygen(users_secret_keys[i]));

        indices.push_back(i + 1);
    }

    std::cout << "Users individually and randomly choose their x and x_dash values.\n";
    std::vector<BICYCL::Mpz> users_secret_x;
    std::vector<BICYCL::Mpz> users_secret_x_dash;

    for (size_t i = 0; i < n; ++i) {
        users_secret_x.push_back(randgen.random_mpz(pp.encrypt_randomness_bound()));
        users_secret_x_dash.push_back(randgen.random_mpz(pp.encrypt_randomness_bound()));
    }

    std::cout << "Users individually compute integer shares of their x and x_dash values.\n";
    std::vector<std::vector<BICYCL::Mpz>> users_x_shares;
    std::vector<std::vector<BICYCL::Mpz>> users_x_dash_shares;

    for (size_t i = 0; i < n; ++i) {
        IntegerPolynomial ip_x = IntegerPolynomial(t, n, pp.secretkey_bound(), users_secret_x[i], randgen);
        IntegerPolynomial ip_x_dash = IntegerPolynomial(t, n, pp.secretkey_bound(), users_secret_x_dash[i], randgen);

        std::vector<BICYCL::Mpz> x_shares;
        std::vector<BICYCL::Mpz> x_dash_shares;

        for (int j = 1; j <= n; ++j) {
            x_shares.push_back(ip_x.evaluate(j));
            x_dash_shares.push_back(ip_x_dash.evaluate(j));
        }

        users_x_shares.push_back(x_shares);
        users_x_dash_shares.push_back(x_dash_shares);
    }

    std::cout << "Setting base h value for compute class group commitments.\n";
    BICYCL::QFI PC_base; 
    pp.power_of_h(PC_base, randgen.random_mpz(pp.encrypt_randomness_bound()));

    ClassGroupCommitmentBuilder cgc(PC_base, pp);

    std::cout << "Users individually compute Gen(pp) data for each user.\n";
    std::vector< // this list stores data of all n users
                std::vector< // this list stores data generated by i-th user for all n users.
                            std::tuple< // this tuple stores data generated by i-th user for j-th user.
                                       BICYCL::QFI, // class group commitment
                                       std::vector<BICYCL::QFI>, // list of c_0 values of encrypted q_ary of share of x.
                                       std::vector<BICYCL::QFI>, // list of c_1 values of encrypted q_ary of share of x.
                                       BICYCL::QFI, // c_0 value of encryption of group element (lifted delta times share).
                                       BICYCL::QFI, // c_1 value of encryption of group element (lifted delta times share).
                                       BIntZKProof // Zero knowledge proof value for round 1 data.
                                       >
                            >
                > users_1st_round_data;

    for (size_t i = 0; i < n; ++i) {
        std::vector<std::tuple<BICYCL::QFI, std::vector<BICYCL::QFI>, std::vector<BICYCL::QFI>, BICYCL::QFI, BICYCL::QFI, BIntZKProof>> user_data;

        for (size_t j = 0; j < n; ++j) {
            // Each i-th user commit their share of x for j-th user

            BICYCL::QFI PC;
            cgc.compute_commitment(PC, users_x_shares[i][j], users_x_dash_shares[i][j]);

            // Each i-th user compute lifted share * delta and encrypts it for j-th user
            BICYCL::Mpz delta_times_share = users_x_shares[i][j];
            BICYCL::Mpz::mul(delta_times_share, delta_times_share, delta);

            BICYCL::QFI gq_to_pow_of_delta_times_share;
            pp.power_of_h(gq_to_pow_of_delta_times_share, delta_times_share);

            BICYCL::QFI c_0, c_1;
            BICYCL::Mpz r = randgen.random_mpz(pp.encrypt_randomness_bound());
            group_element_encryption(c_0, c_1, users_public_keys[j], 
                                     gq_to_pow_of_delta_times_share, 
                                     r, pp);
            
            // Each i-th user decompose their share for j-th user into q-ary and encrypt all the chunks.
            std::vector<BICYCL::Mpz> x_share_q_ary = to_q_ary(users_x_shares[i][j], pp.q());

            std::vector<BICYCL::QFI> c_0_values;
            std::vector<BICYCL::QFI> c_1_values;
            std::vector<BICYCL::Mpz> r_values;

            for (size_t k = 0; k < x_share_q_ary.size(); ++k) {
                r_values.push_back(randgen.random_mpz(pp.encrypt_randomness_bound()));
                BICYCL::CL_HSMqk::CipherText ct = pp.encrypt(users_public_keys[j], 
                                                             BICYCL::CL_HSMqk::ClearText(pp, x_share_q_ary[k]), 
                                                             r_values[k]);

                c_0_values.push_back(ct.c1());
                c_1_values.push_back(ct.c2());
            }

            // Each i-th user computes ZKP for round 1 data compute for j-th user.
            BIntZKProofGenerationInput input {x_share_q_ary, r_values, r, users_x_dash_shares[i][j]};
            BIntZKProofToProve to_prove {users_public_keys[j], PC, c_0_values, c_1_values, c_0, c_1};

            BIntZKProof proof = generate_BInt_proof(input, to_prove, delta, pp, cgc, randgen);

            user_data.push_back(std::make_tuple(PC, c_0_values, c_1_values, c_0, c_1, proof));
        }

        users_1st_round_data.push_back(user_data);
    }

    std::cout << "Users verify each other's data in GenVf().\n";
    IntegerDualPolynomial idp = IntegerDualPolynomial(t, n, pp.secretkey_bound(), randgen);

    std::vector<BICYCL::Mpz> duals;
    for (size_t i = 0; i < n; ++i) {
        duals.push_back(idp.evaluate(i+1, indices));
    }

    for (size_t i = 0; i < n; ++i) {
        std::vector<BICYCL::QFI> commitments_created_by_ith_user;
        for (size_t j = 0; j < n; ++j) {
            commitments_created_by_ith_user.push_back(std::get<0>(users_1st_round_data[i][j]));
        }

        if (IntegerDualPolynomial::verification_in_commitment(commitments_created_by_ith_user, duals, pp) == false) {
            // remove i-th user from list of participants.
            std::cout << "Removed " << i + 1 << "-th user because they generated inconsistent shares.\n";
            auto it = std::find(indices.begin(), indices.end(), i + 1); 
            if (it != indices.end()) { 
                indices.erase(it); 
            } 
        } else {
            for (size_t j = 0; j < n; ++j) {
                BIntZKProofToProve to_prove {users_public_keys[j], 
                                             std::get<0>(users_1st_round_data[i][j]),
                                             std::get<1>(users_1st_round_data[i][j]),
                                             std::get<2>(users_1st_round_data[i][j]),
                                             std::get<3>(users_1st_round_data[i][j]),
                                             std::get<4>(users_1st_round_data[i][j])};

                bool result = verify_BInt_proof(std::get<5>(users_1st_round_data[i][j]), to_prove, delta, pp, cgc);
                if (result == false) {
                    // remove i-th user from list of participants.
                    std::cout << "Removed " << i + 1 << "-th user because their ZKPs didn't verify.\n";
                    auto it = std::find(indices.begin(), indices.end(), i + 1); 
                    if (it != indices.end()) { 
                        indices.erase(it); 
                    } 
                }    
            }
        }
    }

    std::cout << "Each user adds their shares from other users and reveal public component in Reveal().\n";
    std::vector<BICYCL::Mpz> decryption_key_shares;

    std::vector< // this list stores round 2 data of all n users
                std::tuple< // each user has this data
                           BICYCL::QFI, // Xi_cross
                           GDecCLZKProof // GDec-CL ZKP
                          >
               > users_2nd_round_data;

    for (size_t i = 0; i < n; ++i) {
        BICYCL::Mpz x_i(0UL);
        BICYCL::QFI c0_x_cross_i;
        BICYCL::QFI c1_x_cross_i;

        bool first_time = true;

        for (auto j : indices) {
            std::vector<BICYCL::Mpz> decrypted_q_ary;
            std::vector<BICYCL::QFI> c0_values = std::get<1>(users_1st_round_data[j-1][i]);
            std::vector<BICYCL::QFI> c1_values = std::get<2>(users_1st_round_data[j-1][i]);

            for (size_t k = 0; k < c0_values.size(); ++k) {
                BICYCL::QFI fm;

                pp.Cl_G().nupow (fm, c0_values[k], users_secret_keys[i]);
                pp.Cl_Delta().nucompinv (fm, c1_values[k], fm); 

                decrypted_q_ary.push_back(pp.dlog_in_F(fm));
            }

            if (decrypted_q_ary.size() != c0_values.size()) {
                throw std::invalid_argument("Decryption of q_ary failed.\n");
            }

            BICYCL::Mpz x_ji = from_q_ary(decrypted_q_ary, pp.q());
            BICYCL::Mpz::add(x_i, x_i, x_ji);

            if (first_time) {
                c0_x_cross_i = std::get<3>(users_1st_round_data[j-1][i]);
                c1_x_cross_i = std::get<4>(users_1st_round_data[j-1][i]);
            } else {
                BICYCL::QFI c0, c1;
                BICYCL::Mpz r = randgen.random_mpz(pp.encrypt_randomness_bound());

                pp.power_of_h(c0, r);
                users_public_keys[i].exponentiation(pp, c1, r);

                pp.Cl_G().nucomp(c0, c0, std::get<3>(users_1st_round_data[j-1][i]));
                pp.Cl_G().nucomp(c0, c0, c0_x_cross_i);

                pp.Cl_Delta().nucomp(c1, c1, std::get<4>(users_1st_round_data[j-1][i]));
                pp.Cl_Delta().nucomp(c1, c1, c1_x_cross_i);

                c0_x_cross_i = c0;
                c1_x_cross_i = c1;
            }

            first_time = false;
        }
    
        decryption_key_shares.push_back(x_i);
        BICYCL::Mpz delta_times_x_i = delta;
        BICYCL::Mpz::mul(delta_times_x_i, delta_times_x_i, x_i);

        BICYCL::QFI Xi_cross;
        pp.power_of_h(Xi_cross, delta_times_x_i);

        GDecCLZKProofGenerationInput input = {users_secret_keys[i], x_i};
        GDecCLZKProofToProve to_prove = {Xi_cross, c0_x_cross_i, c1_x_cross_i, users_public_keys[i]};

        GDecCLZKProof proof = generate_GDecCL_proof(input, to_prove, delta, pp, randgen);

        users_2nd_round_data.push_back(std::make_tuple(Xi_cross, proof));
    }

    return 0;
}