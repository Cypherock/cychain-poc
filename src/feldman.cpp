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

class PrimeSeq {
public:
    PrimeSeq() : current_prime_index(0) {
        // Precompute a list of prime numbers using the sieve of Eratosthenes up to a reasonable bound
        generate_primes(10000); // Can adjust the bound as needed
    }

    // Resets the sequence, starting again at the first prime
    void reset(long start = 2) {
        current_prime_index = 0;
        while (current_prime_index < primes.size() && primes[current_prime_index] < start) {
            ++current_prime_index;
        }
    }

    // Returns the next prime in the sequence, or 0 if we reach the end of the list
    long next() {
        if (current_prime_index < primes.size()) {
            return primes[current_prime_index++];
        }
        return 0; // No more primes in the sequence
    }

private:
    std::vector<long> primes;
    size_t current_prime_index;

    // Generates primes using the Sieve of Eratosthenes up to the given bound
    void generate_primes(long bound) {
        std::vector<bool> is_prime(bound + 1, true);
        is_prime[0] = is_prime[1] = false;

        for (long p = 2; p * p <= bound; ++p) {
            if (is_prime[p]) {
                for (long i = p * p; i <= bound; i += p) {
                    is_prime[i] = false;
                }
            }
        }

        // Store all the primes in the vector
        for (long p = 2; p <= bound; ++p) {
            if (is_prime[p]) {
                primes.push_back(p);
            }
        }
    }
};

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


class FeldmanVSS_BICYCL {
    static BICYCL::Mpz q;
    static BICYCL::Mpz g;
public:
    vector<BICYCL::Mpz> shares;
    vector<BICYCL::Mpz> commits;
    vector<uchar> cipher;

    FeldmanVSS_BICYCL();

    void print();

    static void init();
    static void load(const BICYCL::Mpz &q, const BICYCL::Mpz &g);

    static BICYCL::Mpz get_q();
    static BICYCL::Mpz get_g();

    static BICYCL::Mpz commit(const BICYCL::Mpz &a);
    static bool verify(const vector<BICYCL::Mpz> commits, int xval, const BICYCL::Mpz &share);

    // Currently relies on the fact that data is at least 16 bytes
    static FeldmanVSS_BICYCL split(const vector<uchar> &data);
    static vector<uchar> reconstruct(const vector<int> &xvec, const vector<BICYCL::Mpz> &shares, const vector<uchar> &cipher);
};

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

// class BICYCL_px which implements polynomial arithmetic using BICYCL library
class BICYCL_px {
  private:
    std::vector<Mpz> coeffs;

  public:
    BICYCL_px() = default;
    explicit BICYCL_px(const std::vector<Mpz>& c) : coeffs(c) {}

    // Degree of the polynomial
    size_t degree() const {
        return coeffs.empty() ? 0 : coeffs.size() - 1;
    }

    // Access coefficient
    Mpz& operator[](size_t i) {
        return coeffs[i];
    }

    const Mpz& operator[](size_t i) const {
        return coeffs[i];
    }

    // Addition of polynomials
    BICYCL_px operator+(const BICYCL_px& other) const {
        size_t max_deg = std::max(degree(), other.degree());
        std::vector<Mpz> result_coeffs(max_deg + 1);

        for (size_t i = 0; i <= max_deg; ++i) {
            if (i <= degree()) result_coeffs[i] = coeffs[i];
            if (i <= other.degree()) BICYCL::Mpz::add(result_coeffs[i], result_coeffs[i], other.coeffs[i]);
        }

        return BICYCL_px(result_coeffs);
    }

    // Subtraction of polynomials
    BICYCL_px operator-(const BICYCL_px& other) const {
        size_t max_deg = std::max(degree(), other.degree());
        std::vector<Mpz> result_coeffs(max_deg + 1);

        for (size_t i = 0; i <= max_deg; ++i) {
            if (i <= degree()) result_coeffs[i] = coeffs[i];
            if (i <= other.degree()) BICYCL::Mpz::sub(result_coeffs[i], result_coeffs[i], other.coeffs[i]);
        }

        return BICYCL_px(result_coeffs);
    }

    // Multiplication of polynomials
    BICYCL_px operator*(const BICYCL_px& other) const {
        size_t result_deg = degree() + other.degree();
        std::vector<Mpz> result_coeffs(result_deg + 1);

        for (size_t i = 0; i <= degree(); ++i) {
            for (size_t j = 0; j <= other.degree(); ++j) {
                BICYCL::Mpz product;
                BICYCL::Mpz::mul(product, coeffs[i], other.coeffs[j]);
                BICYCL::Mpz::add(result_coeffs[i + j], result_coeffs[i + j], product);
            }
        }

        return BICYCL_px(result_coeffs);
    }

    // Evaluate polynomial at a point
    BICYCL::Mpz eval(const BICYCL::Mpz& x) const {
        BICYCL::Mpz result(0UL);
        BICYCL::Mpz x_pow(1UL);

        for (size_t i = 0; i <= degree(); ++i) {
            BICYCL::Mpz term;
            BICYCL::Mpz::mul(term, coeffs[i], x_pow);
            BICYCL::Mpz::add(result, result, term);
            BICYCL::Mpz::mul(x_pow, x_pow, x);
        }

        return result;
    }

    // Print polynomial
    friend std::ostream& operator<<(std::ostream& os, const BICYCL_px& poly) {
        for (size_t i = 0; i <= poly.degree(); ++i) {
            os << poly.coeffs[i] << "x^" << i;
            if (i != poly.degree()) os << " + ";
        }
        return os;
    }

};



void FeldmanVSS_BICYCL::init() {
    // Initialize q and g with some default values

    // q = BICYCL::Mpz("some_large_prime_value");
    // g = BICYCL::Mpz("some_generator_value");

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
