#include "CryptoUtils.h"

extern "C"
{
#include "bignum.h"
}

unsigned char CryptoUtils::hex_char_to_int(char c) {
    if (c >= '0' && c <= '9')
        return c - '0';
    if (c >= 'a' && c <= 'f')
        return 10 + c - 'a';
    if (c >= 'A' && c <= 'F')
        return 10 + c - 'A';
    return 0;
}

void CryptoUtils::hex_string_to_byte_array(const char *hex_string, uint8_t *byte_array, size_t byte_array_size) {
    size_t hex_string_length = strlen(hex_string);
    if (hex_string_length % 2 != 0) {
        printf("Hex string must have an even number of characters.\n");
        return;
    }
    if (byte_array_size < hex_string_length / 2) {
        printf("Byte array is too small to hold the converted hex string.\n");
        return;
    }
    for (size_t i = 0; i < hex_string_length; i += 2) {
        byte_array[i / 2] = (hex_char_to_int(hex_string[i]) << 4) + hex_char_to_int(hex_string[i + 1]);
    }
}

void CryptoUtils::byte_to_hex(unsigned char byte, char hex[3]) {
    const char hex_chars[] = "0123456789abcdef";
    hex[0] = hex_chars[(byte >> 4) & 0x0F];
    hex[1] = hex_chars[byte & 0x0F];
    hex[2] = '\0';
}

void CryptoUtils::byte_array_to_hex_string(const uint8_t *byte_array, size_t byte_array_size, char *hex_string) {
    hex_string[0] = '\0'; // Start with an empty string
    for (size_t i = 0; i < byte_array_size; i++) {
        char hex[3];
        byte_to_hex(byte_array[i], hex);
        strcat(hex_string, hex);
    }
}

void CryptoUtils::print_bytes(uint8_t *bytes, size_t len) {
    for (size_t i = 0; i < len; i++) {
        printf("%02x", bytes[i]);
    }
    printf("\n");
}

void CryptoUtils::print_bignum256(bignum256 *bn) {
    uint8_t bytes[32];
    bn_write_be(bn, bytes);
    for (int i = 0; i < 32; i++) {
        printf("%02x", bytes[i]);
    }
    printf("\n");
}

void CryptoUtils::rand_bytes(uint8_t *bytes, size_t len) {
    if (RAND_bytes(bytes, len) != 1) {
        fprintf(stderr, "Error generating random bytes\n");
        exit(EXIT_FAILURE);
    }
}

void CryptoUtils::rand_bignum256(const bignum256 *curve_order, bignum256 *res) {
    uint8_t random_bytes[32];
    rand_bytes(random_bytes, 32);
    bn_read_be(random_bytes, res);
    while (!bn_is_less(res, curve_order)) {
        bn_mod(res, curve_order);
    }
}
