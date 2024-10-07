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
using namespace BICYCL;

BICYCL::Mpz mod_exp(BICYCL::Mpz base, BICYCL::Mpz exp, BICYCL::Mpz mod) {
    BICYCL::Mpz result(1UL);
    BICYCL::Mpz::mod(base, base, mod);
    while (BICYCL::Mpz::operator>(exp, BICYCL::Mpz(0UL))) {
        BICYCL::Mpz temp_exp;
        BICYCL::Mpz::mod(temp_exp, exp, BICYCL::Mpz(2UL));
        if(BICYCL::Mpz::operator==(temp_exp, BICYCL::Mpz(1UL))) {
            BICYCL::Mpz::mul(result, result, base);
            BICYCL::Mpz::mod(result, result, mod);

        }
        BICYCL::Mpz::divby2(exp, exp);
        BICYCL::Mpz::mul(base, base, base);
        BICYCL::Mpz::mod(base, base, mod);
    }
    return result;
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
    for (int i = 1; i <= n; ++i) {
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
