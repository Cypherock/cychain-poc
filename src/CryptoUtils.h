#pragma once

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <openssl/rand.h>

extern "C"
{
#include "bignum.h"
}

class CryptoUtils {
public:
    static unsigned char hex_char_to_int(char c);
    static void hex_string_to_byte_array(const char *hex_string, uint8_t *byte_array, size_t byte_array_size);
    static void byte_to_hex(unsigned char byte, char hex[3]);
    static void byte_array_to_hex_string(const uint8_t *byte_array, size_t byte_array_size, char *hex_string);
    static void print_bytes(uint8_t *bytes, size_t len);
    static void print_bignum256(bignum256 *bn);
    static void rand_bytes(uint8_t *bytes, size_t len);
    static void rand_bignum256(const bignum256 *curve_order, bignum256 *res);
};
