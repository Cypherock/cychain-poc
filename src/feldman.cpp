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
#include <vector>
#include <future>

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

    // define a var of type PrimeSeq using BICYCL library & trezor-crypto library
    

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
            if (r == p - r - 1) {
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



// Generates random coefficients for the secret sharing polynomial
std::vector<BICYCL::Mpz> generate_coefficients(int k, BICYCL::Mpz secret, BICYCL::Mpz q) {
    std::vector<BICYCL::Mpz> coefficients(k);
    coefficients[0] = secret;  // First coefficient is the secret
    BICYCL::Mpz::mod(coefficients[0], coefficients[0], q);
    BICYCL::Mpz::mod(coefficients[0], coefficients[0], q);
    BICYCL::RandGen randgen = BICYCL::RandGen(q);
    for (int i = 1; i < k; ++i) {
        coefficients[i] = randgen.random_mpz(q);  // Random coefficients
    }
    return coefficients;
}

// Evaluate polynomial at x for BICYCL::Mpz
BICYCL::Mpz evaluate_polynomial(const std::vector<BICYCL::Mpz>& coefficients, BICYCL::Mpz x, BICYCL::Mpz q) {
    BICYCL::Mpz result(0UL);
    BICYCL::Mpz x_pow(1UL);
    for (BICYCL::Mpz coeff : coefficients) {
        BICYCL::Mpz temp;
        BICYCL::Mpz::mul(temp, coeff, x_pow);
        BICYCL::Mpz::add(result, result, temp);
        BICYCL::Mpz::mod(result, result, q);
        BICYCL::Mpz::mul(x_pow, x_pow, x);
        BICYCL::Mpz::mod(x_pow, x_pow, q);
    }
    return result;
}

// Generate shares and verification values
void generate_shares_and_verifications(int n, int k, BICYCL::Mpz secret, BICYCL::Mpz q, BICYCL::Mpz g,
                                       std::vector<std::pair<BICYCL::Mpz, BICYCL::Mpz>>& shares,
                                       std::vector<BICYCL::Mpz>& verification_values) {
    // Generate the polynomial coefficients
    std::vector<BICYCL::Mpz> coefficients = generate_coefficients(k, secret, q);

    // Generate verification values A_j = g^a_j mod q
    verification_values.resize(k);
    for (int j = 0; j < k; ++j) {
        verification_values[j] = mod_exp(g, coefficients[j], q);
    }

    // Generate shares (x_i, f(x_i)) for each player P_i
    for (long i = 1; i <= n; ++i) {
        BICYCL::Mpz x_i(i);
        BICYCL::Mpz y_i = evaluate_polynomial(coefficients, x_i, q);
        shares.push_back({x_i, y_i});
    }
}

// Verify a share using the public verification values for BICYCL::Mpz
bool verify_share(const std::pair<BICYCL::Mpz, BICYCL::Mpz>& share, const std::vector<BICYCL::Mpz>& verification_values, BICYCL::Mpz q, BICYCL::Mpz g) {
    BICYCL::Mpz x_i = share.first;
    BICYCL::Mpz f_x_i = share.second;

    // Compute V_i = product(A_j ^ x_i^j) mod q
    BICYCL::Mpz V_i(1UL);
    BICYCL::Mpz x_pow(1UL);

    for (int j = 0; j < verification_values.size(); ++j) {
        BICYCL::Mpz temp;
        BICYCL::Mpz::pow_mod(temp, verification_values[j], x_pow, q);
        BICYCL::Mpz::mul(V_i, V_i, temp);
        BICYCL::Mpz::mod(V_i, V_i, q);
        BICYCL::Mpz::mul(x_pow, x_pow, x_i);
        BICYCL::Mpz::mod(x_pow, x_pow, q);
    }

    // Compute g^f(x_i) mod q
    BICYCL::Mpz V_i_prime = mod_exp(g, f_x_i, q);

    return V_i == V_i_prime;
}

// Main function to run the Feldman's VSS scheme
int feldman() {
    BICYCL::Mpz secret(1234UL);  // Secret to be shared
    int n = 5;                   // Number of players
    int k = 3;                   // Threshold number of shares
    BICYCL::Mpz q(7919UL);       // A large prime modulus
    BICYCL::Mpz g(2UL);          // Generator for Z_q*

    // Generate shares and verification values
    std::vector<std::pair<BICYCL::Mpz, BICYCL::Mpz>> shares;
    std::vector<BICYCL::Mpz> verification_values;
    generate_shares_and_verifications(n, k, secret, q, g, shares, verification_values);

    // Output the shares and verification values
    std::cout << "Shares (x_i, f(x_i)):" << std::endl;
    for (const auto& share : shares) {
        std::cout << "(" << share.first << ", " << share.second << ")" << std::endl;
    }

    std::cout << "Verification values A_j:" << std::endl;
    for (const BICYCL::Mpz& value : verification_values) {
        std::cout << value << " ";
    }
    std::cout << std::endl;

    // Verify each share
    for (const auto& share : shares) {
        bool valid = verify_share(share, verification_values, q, g);
        std::cout << "Share (" << share.first << ", " << share.second << ") is "
                  << (valid ? "valid" : "invalid") << std::endl;
    }

    return 0;
}

