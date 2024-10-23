// feldman.h

#ifndef FELDMAN_H
#define FELDMAN_H

#include <string>
#include <vector>
#include <map>
#include <random>
#include <openssl/rand.h>
#include <stdint.h>
#include <future>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <functional>
#include <atomic>

#include "ThreadPool.h"
#include "CryptoUtils.h"
#include "bicycl.hpp"

extern "C" {
#include "ecdsa.h"
#include "bignum.h"
#include "curves.h"
#include "bip32.h"
#include "sha2.h"
}

#define uchar unsigned char

using std::string;
using namespace std;
using namespace BICYCL;


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

    bool prove_secret_correctness(const vector<BICYCL::Mpz> &commits, const vector<int> &xvec, const vector<BICYCL::Mpz> &sharesl);
    bool verify_secret_correctness(const vector<BICYCL::Mpz> &commits, const vector<int> &xvec, const vector<BICYCL::Mpz> &shares);
    bool prove_share_correctness(const BICYCL::Mpz &commit, const BICYCL::Mpz &share, int xval);
    bool verify_share_correctness(const BICYCL::Mpz &commit, const BICYCL::Mpz &share, int xval);
    bool prove_same_secret(const vector<uchar> &original_secret, const vector<uchar> &reconstructed_secret); 
    bool verify_same_secret(const vector<uchar> &original_secret, const vector<uchar> &reconstructed_secret);
};
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

#endif // FELDMAN_H