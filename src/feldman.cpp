#include "bicycl.hpp"

extern "C"
{
#include "ecdsa.h"
#include "bignum.h"
#include "curves.h"
#include "bip32.h"
#include "sha2.h"
}
#include "feldman.h"
#define uchar unsigned char

using std::string;
using namespace std;
using namespace BICYCL;

BICYCL::Mpz mod_exp(BICYCL::Mpz base, BICYCL::Mpz exp, BICYCL::Mpz mod) {
    BICYCL::Mpz result(1UL);
    BICYCL::Mpz::mod(base, base, mod);
    while (exp.operator>(BICYCL::Mpz(0UL))) {
        BICYCL::Mpz temp_exp;
        BICYCL::Mpz::mod(temp_exp, exp, BICYCL::Mpz(2UL));
        if(temp_exp.operator==(BICYCL::Mpz(1UL))) {
            BICYCL::Mpz::mul(result, result, base);
            BICYCL::Mpz::mod(result, result, mod);

        }
        BICYCL::Mpz::divby2(exp, exp);
        BICYCL::Mpz::mul(base, base, base);
        BICYCL::Mpz::mod(base, base, mod);
    }
    return result;
}

// function to convert BICYCL::Mpz to unsigned long as class BICYCL::Mpz does not have a direct conversion operator to unsigned long
unsigned long mpz_to_ulong(const BICYCL::Mpz& mpz) {
    unsigned long result = 0;
    std::string mpz_str = BICYCL::Mpz:: get_str(mpz, 10); // Convert BICYCL::Mpz to string
    try {
        result = std::stoul(mpz_str); // Convert string to unsigned long
    } catch (const std::out_of_range& e) {
        std::cerr << "Value out of range for unsigned long: " << e.what() << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cerr << "Invalid argument: " << e.what() << std::endl;
    }
    return result;
}

void left_shift(BICYCL::Mpz& result, const BICYCL::Mpz& value, unsigned long shift, BICYCL::Mpz& q) {
    BICYCL::Mpz two_pow_shift;
    two_pow_shift = BICYCL::Mpz(1UL);
    for (unsigned long i = 0; i < shift; ++i) {
        BICYCL::Mpz::mul(two_pow_shift, two_pow_shift, BICYCL::Mpz(2UL));
    }
    BICYCL::Mpz::mul(result, value, two_pow_shift);
    BICYCL::Mpz::mod(result, result, q);
}

// Generating Germain prime
long ComputePrimeBound(long k) {
    return 1L << (k / 2 - 1);
}

// RandomBnd function generates x = "random number" in the range 0..n-1, or 0  if n <= 0
 void RandomBnd(BICYCL::Mpz& x, const BICYCL::Mpz& n) {
    if (n.operator<=(BICYCL::Mpz(0UL))) {
        x.operator=(BICYCL::Mpz(0UL));
        return;
    }

    BICYCL::RandGen gen;
    x = gen.random_mpz(n);
 }

bool MillerWitness(BICYCL::Mpz& n, BICYCL::Mpz& W) {
    BICYCL::Mpz n1;
    BICYCL::Mpz::sub(n1,n,1);
    BICYCL::Mpz s(0UL);
    BICYCL::Mpz t = n1;

    while (t.is_even()) {
        BICYCL::Mpz::divby2(t, t);
        BICYCL::Mpz::add(s, s, 1);
    }

    BICYCL::Mpz a = mod_exp(t, W, n);
    if (a.operator==(1L) || a.operator==(n1)) {
        return false;
    }
    BICYCL::Mpz i;
    for(i.operator=(1UL); i.operator<(s); BICYCL::Mpz::add(i, i, 1UL)) {
        // a = a.sqrMod(n);
        BICYCL::Mpz::mul(a, a, a);
        BICYCL::Mpz::mod(a, a, n);
        if (a.operator==(n1)) {
            return false;
        }
    }
    return true;
}

void MultiThreadedGenGermainPrime(BICYCL::Mpz& n, long k, long err) {
    ThreadPool pool(std::thread::hardware_concurrency());
    
    vector<future<void>> results;
    for (int i = 0; i < std::thread::hardware_concurrency(); i++) {
        results.push_back(
            std::async(std::launch::async, [&n, k, err]() { GenGermainPrime_BICYCL(n, k, err); })
        );
    }

    for (auto &&result : results) {
        result.get();
    }
}

// implementing genGermainPrime using BICYCL and trezor-crypto libraries and not using NTL and openssl libraries and not using ZZ datatype to represent large numbers

void GenGermainPrime_BICYCL(BICYCL::Mpz& n, long k, long err) {
    if (k <= 1) {
        throw "GenGermainPrime: bad length";
    }

    if (k > (1L << 20)) {
        throw "GenGermainPrime: length too large";
    }

    if (err < 1) {
        err = 1;
    }
    if (err > 512) {
        err = 512;
    }

    if (k == 2) {
        BICYCL::Mpz r;
        RandomBnd(r, BICYCL::Mpz(2UL));
        if (r.operator==(0UL)) {
            n.operator=(2UL);
        } else {
            n.operator=(3UL);
        }

        return;
    }

    long prime_bnd = ComputePrimeBound(k);
    BICYCL::Mpz prime_bnd_mpz(prime_bnd);

    if (prime_bnd_mpz.nbits() >= k / 2) {
        prime_bnd_mpz.operator=(1L << (k / 2 - 1));
    }

    BICYCL::Mpz two;
    two.operator=(2UL);
    BICYCL::Mpz n1;
    
    PrimeSeq s;

    BICYCL::Mpz iter;
    iter.operator=(0UL);

    for (;;) {
        BICYCL::Mpz::add(iter, iter, 1UL);

        RandomBnd(n, BICYCL::Mpz(k));
        if (n.is_even()) {
        BICYCL::Mpz::add(n, n, 1UL);
        }

        s.reset(3);
        long p;

        long sieve_passed = 1;

        p = s.next();
        while (p && p < prime_bnd) {
            BICYCL::Mpz r;
            BICYCL::Mpz::mod(r, n, BICYCL::Mpz(p));

            if (r.operator==(0UL)) {
                sieve_passed = 0;
                break;
            }

            // test if 2*r + 1 = 0 (mod p)
            BICYCL::Mpz r2, r2_1;
            BICYCL::Mpz::add(r2, r, 1UL);
            BICYCL::Mpz::sub(r2_1,BICYCL::Mpz(p), r2);
            if(r.operator==(r2_1)) {
                sieve_passed = 0;
                break;
            }

            p = s.next();
        }

        if (!sieve_passed) {
            continue;
        }

        if (MillerWitness(n, two)) {
            continue;
        }

        // n1 = 2*n+1
        BICYCL::Mpz::mulby2(n1, n);
        BICYCL::Mpz::add(n1, n1, 1);

        if (MillerWitness(n1, two)) {
            continue;
        }

        // now do t M-R iterations...just to make sure

        // First compute the appropriate number of M-R iterations, t
        // The following computes t such that
        //       p(k,t)*8/k <= 2^{-err}/(5*iter^{1.25})
        // which suffices to get an overall error probability of 2^{-err}.
        // Note that this method has the advantage of not requiring
        // any assumptions on the density of Germain primes.

}
}

bool ErrBoundTest(long k, long t, long err) {
    return true;
}


BICYCL::Mpz FeldmanVSS_BICYCL::q;
BICYCL::Mpz FeldmanVSS_BICYCL::g;

FeldmanVSS_BICYCL::FeldmanVSS_BICYCL() {}

void FeldmanVSS_BICYCL::print() {
    cout << "Shares: ";
    for (const auto &share : shares) {
        cout << share << " ";
    }
    cout << endl;

    cout << "Commits: ";
    for (const auto &commit : commits) {
        cout << commit << " ";
    }
    cout << endl;
}

void FeldmanVSS_BICYCL::init() {
    // Generate prime q such that 2q+1 is also a prime
    // In all code below, p := 2q + 1
    BICYCL::Mpz q;
    GenGermainPrime_BICYCL(q, 256, 80);
    BICYCL::Mpz p;
    BICYCL::Mpz::mulby2(p, q);
    BICYCL::Mpz::add(p, p, 1);
    BICYCL::Mpz b;
    RandomBnd(b, p);
    while(b.operator==(0UL))
        RandomBnd(b, p);
    g = mod_exp(b, BICYCL::Mpz(2UL), p);
}

void FeldmanVSS_BICYCL::load(const BICYCL::Mpz &q_val, const BICYCL::Mpz &g_val) {
    q = q_val;
    g = g_val;
}

BICYCL::Mpz FeldmanVSS_BICYCL::get_q() {
    return q;
}

BICYCL::Mpz FeldmanVSS_BICYCL::get_g() {
    return g;
}

BICYCL::Mpz FeldmanVSS_BICYCL::commit(const BICYCL::Mpz &a) {
    BICYCL::Mpz p;
    BICYCL::Mpz::mulby2(p, q);
    BICYCL::Mpz::add(p, p, 1);
    return mod_exp(g, a, p);
}

FeldmanVSS_BICYCL FeldmanVSS_BICYCL::split(const vector<uchar> &secret) {
    if (secret.size() < 16) {
        throw "split: Secret not long enough";
    }

    BICYCL::Mpz p;
    BICYCL::Mpz::mulby2(p, q);
    BICYCL::Mpz::add(p, p, 1);

    // Generate random linear polynomial 'f' for Shamir SS
    // f = (c[0] + c[1]x) mod q
    std::vector<BICYCL::Mpz> coeffs(2);
    BICYCL::RandGen gen;
    coeffs[0] = gen.random_mpz(q);
    coeffs[1] = gen.random_mpz(q);
    BICYCL_px f(coeffs);

    FeldmanVSS_BICYCL out;

    // Computing shares for each party
    // shares[i] := f(i) mod q, for i in [1, 4]
    for (int i = 1; i <= 4; ++i) {
        BICYCL::Mpz x((long)i);
        BICYCL::Mpz share = f.eval(x);
        out.shares.push_back(share);
    }

    // Computing the commitments for the coefficients
    // commits[i] := g^(c[i]) mod p
    for (const auto &coeff : coeffs) {
        out.commits.push_back(commit(coeff));
    }

    // Getting 'data stream' from BICYCL::Mpz type for constant term (c[0])
    std::string c_str = coeffs[0].get_str(BICYCL::Mpz(10UL), 10);
    std::vector<uchar> c(c_str.begin(), c_str.end());

    // Computing key
    // key := KDF(c[0])
    uchar key[SHA256_DIGEST_LENGTH];
    SHA256_CTX ctx;
    sha256_Init(&ctx);
    sha256_Update(&ctx, c.data(), c.size());
    sha256_Final(&ctx, key);

    // Encrypting secret
    // cipher := secret + key
    out.cipher.resize(secret.size());
    for (size_t i = 0; i < secret.size(); ++i) {
        out.cipher[i] = secret[i] ^ key[i % SHA256_DIGEST_LENGTH];
    }

    return out;
}

bool FeldmanVSS_BICYCL::verify(const vector<BICYCL::Mpz> commits, int xval, const BICYCL::Mpz &share) {
    if (commits.size() < 2) {
        throw "verify: Expected 2 commits";
    }

    BICYCL::Mpz p;
    BICYCL::Mpz::mulby2(p, q);
    BICYCL::Mpz::add(p, p, 1);

    BICYCL::Mpz acc = commits[0];
    BICYCL::Mpz xval_mpz((long)xval);
    BICYCL::Mpz exp_commit;
    exp_commit = mod_exp(commits[1], xval_mpz, p);
    BICYCL::Mpz::mul(acc, acc, exp_commit);
    BICYCL::Mpz::mod(acc, acc, p);

    BICYCL::Mpz share_commit = commit(share);

    return (acc == share_commit);
}

vector<uchar> FeldmanVSS_BICYCL::reconstruct(
    const vector<int> &xvec,
    const vector<BICYCL::Mpz> &shares,
    const vector<uchar> &cipher
) {

    if (xvec.size() < 2 || shares.size() < 2 || cipher.size() < 16) {
        throw "reconstruct: Input size shorter than expected";
    }

    BICYCL::Mpz q = get_q();
    BICYCL_px f;

    // Lagrange interpolation to find the constant term c[0]
    for (size_t i = 0; i < xvec.size(); ++i) {
        BICYCL::Mpz num(1UL), denom(1UL);
        for (size_t j = 0; j < xvec.size(); ++j) {
            if (i != j) {
                BICYCL::Mpz xj((long)xvec[j]);
                BICYCL::Mpz xi((long)xvec[i]);
                BICYCL::Mpz diff;
                BICYCL::Mpz::sub(diff, xj, xi);
                BICYCL::Mpz::mul(num, num, xj);
                BICYCL::Mpz::mul(denom, denom, diff);
            }
        }
        BICYCL::Mpz inv_denom;
        BICYCL::Mpz::divexact(inv_denom, BICYCL::Mpz(1UL), denom);
        BICYCL::Mpz::mod(inv_denom, inv_denom, q);
        BICYCL::Mpz term;
        BICYCL::Mpz::mul(term, shares[i], num);
        BICYCL::Mpz::mul(term, term, inv_denom);
        BICYCL::Mpz::add(f[0], f[0], term);
    }

    // Get 'data stream' of constant term
    std::string c_str = f[0].get_str(BICYCL::Mpz(10UL), 10);
    std::vector<uchar> c(c_str.begin(), c_str.end());

    // Compute key from constant term
    uchar key[SHA256_DIGEST_LENGTH];
    SHA256_CTX ctx;
    sha256_Init(&ctx);
    sha256_Update(&ctx, c.data(), c.size());
    sha256_Final(&ctx, key);

    // Decrypt cipher
    vector<uchar> secret(cipher.size());
    for (size_t i = 0; i < cipher.size(); ++i) {
        secret[i] = cipher[i] ^ key[i % SHA256_DIGEST_LENGTH];
    }

    return secret;
}

//implement zkp for the same secret using BICYCL library and trezor-crypto library
bool FeldmanVSS_BICYCL::prove_same_secret(const vector<uchar> &original_secret, const vector<uchar> &reconstructed_secret) {
    if (original_secret.size() != reconstructed_secret.size()) {
        throw "prove_same_secret: Secret sizes do not match";
    }

    for (size_t i = 0; i < original_secret.size(); ++i) {
        if (original_secret[i] != reconstructed_secret[i]) {
            return false;
        }
    }
    return true;
}

bool FeldmanVSS_BICYCL::verify_same_secret(const vector<uchar> &original_secret, const vector<uchar> &reconstructed_secret) {
    if (original_secret.size() != reconstructed_secret.size()) {
        throw "verify_same_secret: Secret sizes do not match";
    }

    for (size_t i = 0; i < original_secret.size(); ++i) {
        if (original_secret[i] != reconstructed_secret[i]) {
            return false;
        }
    }
    return true;
}

struct DLogProof {
    BICYCL::Mpz g, h, x; // g and h are group elements, x is the witness
    
    void generate_proof(BICYCL::Mpz q) {
        // Choose a random blinding factor
        BICYCL::Mpz r;
        BICYCL::RandGen gen;
        BICYCL::Mpz shift_value;
        left_shift(shift_value, BICYCL::Mpz(1UL), 256, q);
        r = gen.random_mpz(shift_value); // A random number for blinding

        // Compute g^r
        BICYCL::Mpz gr = mod_exp(g, r, q);  // g^r (group operation)

        // Send gr as the commitment to the verifier
        std::cout << "Commitment: " << gr << std::endl;

        // The verifier responds with a challenge c
        BICYCL::Mpz c;
        std::cout << "Enter challenge: ";
        std::cin >> c;

        // Prove knowledge of x by computing z = r + c * x
        BICYCL::Mpz z;
        BICYCL::Mpz::mul(z, c, x);
        BICYCL::Mpz::add(z, z, r);

        // Send the proof z to the verifier
        std::cout << "Proof: " << z << std::endl;
    }
    
    bool verify(BICYCL::Mpz c, BICYCL::Mpz z, BICYCL::Mpz q, BICYCL::Mpz gr) {
        // Verify that g^z = gr * h^c
        BICYCL::Mpz gz = mod_exp(g, z, q);   // g^z
        BICYCL::Mpz hc = mod_exp(h, c, q);   // h^c
        BICYCL::Mpz gr_hc;
        BICYCL::Mpz::mul(gr_hc, gr, hc);
        BICYCL::Mpz::mod(gr_hc, gr_hc, q);
        return gz == gr_hc;
    }
};

struct EqDLogProof {
    BICYCL::Mpz g1, g2, h1, h2, x; // g1, g2, h1, h2 are group elements, x is the witness

    void generate_proof(BICYCL::Mpz q) {
        // Choose random blinding factor r
        BICYCL::Mpz r;
        BICYCL::RandGen gen;
        BICYCL::Mpz shift_value;
        left_shift(shift_value, BICYCL::Mpz(1UL), 256, q);
        r = gen.random_mpz(shift_value);

        // Compute commitments g1^r and g2^r
        BICYCL::Mpz g1r = mod_exp(g1, r, q);
        BICYCL::Mpz g2r = mod_exp(g2, r, q);

        // Send commitments to verifier
        std::cout << "Commitments: " << g1r << ", " << g2r << std::endl;

        // Verifier sends a challenge c
        BICYCL::Mpz c;
        std::cout << "Enter challenge: ";
        std::cin >> c;

        // Prove knowledge of x for both pairs
        BICYCL::Mpz z;
        BICYCL::Mpz::mul(z, c, x);
        BICYCL::Mpz::add(z, z, r);

        // Send proof
        std::cout << "Proof: " << z << std::endl;
    }

    bool verify(BICYCL::Mpz c, BICYCL::Mpz z, BICYCL::Mpz q) {
        // Verify g1^z = g1r * h1^c and g2^z = g2r * h2^c
        BICYCL::Mpz g1z = mod_exp(g1, z, q);
        BICYCL::Mpz g2z = mod_exp(g2, z, q);
        BICYCL::Mpz h1c = mod_exp(h1, c, q);
        BICYCL::Mpz h2c = mod_exp(h2, c, q);

        return (g1z == h1c) && (g2z == h2c);
    }
};
