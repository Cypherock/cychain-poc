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
#include "PolicyContract.h"

#include "bicycl.hpp"

extern "C"
{
#include "ecdsa.h"
#include "bignum.h"
#include "curves.h"
#include "bip32.h"
}

using std::string;

class GroupsPolicyContract : public PolicyContract
{
public:
    // Implementing the checkCompliance method
    bool check_compliance(uint32_t threshold, std::string message_hash_hex_str) override
    {
        size_t size = message_hash_hex_str.length() / 2;
        uint8_t message_bytes[size];

        CryptoUtils::hex_string_to_byte_array(message_hash_hex_str.c_str(), message_bytes, size);

        if (size == 0)
        {
            return false;
        }

        // For PoC, we just check if the first byte of message hash is 0x00
        return message_bytes[0] == 0x01 && threshold != 0;
    }
};

class AccessControlContract
{
private:
    BICYCL::CL_HSMqk pp;
    BICYCL::CL_HSMqk::PublicKey ekv;

    // validators and group will broadcast their sig-gen rounds data to mempool
    std::map<std::string, std::tuple<PolicyContract *, BICYCL::CL_HSMqk::PublicKey, BICYCL::CL_HSMqk::CipherText, std::string>> group_registry;
    std::map<std::string, std::tuple<std::string, BICYCL::CL_HSMqk::CipherText, curve_point, std::string>> round1_data;
    std::map<std::string, std::tuple<std::string, BICYCL::CL_HSMqk::CipherText, BICYCL::CL_HSMqk::CipherText, BICYCL::CL_HSMqk::CipherText>> round2_data;
    std::map<std::string, std::tuple<std::string, BICYCL::CL_HSMqk::CipherText>> round3_data;
    std::map<std::string, std::tuple<std::string, BICYCL::CL_HSMqk::ClearText>> round4_data;
    bool is_setup = false;

public:
    AccessControlContract(BICYCL::CL_HSMqk pp) : pp(pp),
                                                 ekv(BICYCL::CL_HSMqk::PublicKey(
                                                     pp, BICYCL::CL_HSMqk::SecretKey(pp, (BICYCL::RandGen()).random_mpz(pp.secretkey_bound())))) {}

    // Method to set the validator network's encryption key
    void set_encryption_key(const BICYCL::CL_HSMqk::PublicKey &key)
    {
        if (is_setup == true)
        {
            throw std::runtime_error("Access Control Contract already setup");
        }

        ekv = key;
        is_setup = true;
    }

    // Method to get the validator network's encryption key
    BICYCL::CL_HSMqk::PublicKey get_encryption_key()
    {
        if (is_setup == false)
        {
            throw std::runtime_error("Access Control Contract not setup");
        }
        return ekv;
    }

    // Method to register a user group with a policy, public key and cipher text
    std::string register_user_group(PolicyContract *policy, const BICYCL::CL_HSMqk::PublicKey &public_key, const BICYCL::CL_HSMqk::CipherText &cipher_text, uint8_t group_pub_key[65])
    {
        if (is_setup == false)
        {
            throw std::runtime_error("Access Control Contract not setup");
        }

        std::cout << "Registering group on the Access Control Contract...\n";
        // Generate a random unique ID for the group
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(100000, 999999);
        std::string group_id = "group_" + std::to_string(dis(gen));

        std::cout << "Group registered with group ID: " << group_id << std::endl;

        char hex_c_str[131] = {0};
        CryptoUtils::byte_array_to_hex_string(group_pub_key, 65, hex_c_str);

        // Store the policy, public key, and cipher text in the map under the generated ID
        group_registry.insert(std::make_pair(group_id, std::make_tuple(policy, public_key, cipher_text, std::string(hex_c_str))));

        // Return the generated ID
        return group_id;
    }

    // Method to retrieve the PolicyContract for a given group ID
    PolicyContract *get_policy(const std::string &group_id)
    {
        if (is_setup == false)
        {
            throw std::runtime_error("Access Control Contract not setup");
        }

        if (group_registry.find(group_id) != group_registry.end())
        {
            return std::get<0>(group_registry.at(group_id));
        }
        else
        {
            return nullptr; // Return nullptr if the group ID is not found
        }
    }

    // Method to retrieve the PublicKey for a given group ID
    BICYCL::CL_HSMqk::PublicKey get_public_key(const std::string &group_id)
    {
        if (is_setup == false)
        {
            throw std::runtime_error("Access Control Contract not setup");
        }

        if (group_registry.find(group_id) != group_registry.end())
        {
            return std::get<1>(group_registry.at(group_id));
        }
        else
        {
            throw std::runtime_error("Group ID not found"); // Throw an error if the group ID is not found
        }
    }

    // Method to retrieve the CipherText for a given group ID
    BICYCL::CL_HSMqk::CipherText get_cipher_text(const std::string &group_id)
    {
        if (is_setup == false)
        {
            throw std::runtime_error("Access Control Contract not setup");
        }

        if (group_registry.find(group_id) != group_registry.end())
        {
            return std::get<2>(group_registry.at(group_id));
        }
        else
        {
            throw std::runtime_error("Group ID not found"); // Throw an error if the group ID is not found
        }
    }

    // Method to retreive the verification key for a given group ID
    std::string get_verification_key(const std::string &group_id)
    {
        if (is_setup == false)
        {
            throw std::runtime_error("Access Control Contract not setup");
        }

        if (group_registry.find(group_id) != group_registry.end())
        {
            return std::get<3>(group_registry.at(group_id));
        }
        else
        {
            throw std::runtime_error("Group ID not found"); // Throw an error if the group ID is not found
        }
    }

    void add_r1_data(std::vector<std::tuple<BICYCL::CL_HSMqk::CipherText, curve_point>> &r1_data,
                     size_t start_index,
                     size_t end_index,
                     BICYCL::CL_HSMqk::PublicKey ektp,
                     std::tuple<BICYCL::CL_HSMqk::CipherText, curve_point> &result) {
        BICYCL::RandGen randgen;
        const ecdsa_curve *curve = get_curve_by_name(SECP256K1_NAME)->params;
        
        // assign first point to kg_dot_G and then add all the others
        curve_point kg_dot_G = std::get<1>(r1_data.at(start_index));

        // assign first ciphertext to enc_kg and then add all the others
        BICYCL::CL_HSMqk::CipherText enc_kg = std::get<0>(r1_data.at(start_index));

        for (size_t i = start_index + 1; i < end_index; ++i)
        {
            point_add(curve, &std::get<1>(r1_data.at(i)), &kg_dot_G);
            enc_kg = pp.add_ciphertexts(ektp, enc_kg, std::get<0>(r1_data.at(i)), randgen);
        }

        result = std::make_tuple(enc_kg, kg_dot_G);
    }

    std::string upload_ug_r1_data(const std::string &group_id, uint32_t threshold,
                                  std::string message_hash_hex_str,
                                  std::vector<std::tuple<BICYCL::CL_HSMqk::CipherText, curve_point>> r1_data)
    {
        if (is_setup == false)
        {
            throw std::runtime_error("Access Control Contract not setup");
        }

        PolicyContract *group_policy = get_policy(group_id);

        std::cout << "ACC checks the threshold and message against group's policy.\n";
        if (group_policy->check_compliance(threshold, message_hash_hex_str) == false)
        {
            throw std::runtime_error("Signature generation request rejected.");
        }

        std::cout << "Policy compliance check successful.\n";
        if (r1_data.size() == 0)
        {
            throw std::runtime_error("No round 1 data to compute on.");
        }

        const ecdsa_curve *curve = get_curve_by_name(SECP256K1_NAME)->params;
        BICYCL::RandGen randgen;
        BICYCL::CL_HSMqk::PublicKey ektp = get_public_key(group_id);

        // assign first point to kg_dot_G and then add all the others
        curve_point kg_dot_G = std::get<1>(r1_data.at(0));

        // assign first ciphertext to enc_kg and then add all the others
        BICYCL::CL_HSMqk::CipherText enc_kg = std::get<0>(r1_data.at(0));

        std::cout << "Combining all round 1 values of user group...\n";

        for (size_t i = 1; i < r1_data.size(); ++i)
        {
            point_add(curve, &std::get<1>(r1_data.at(i)), &kg_dot_G);
            enc_kg = pp.add_ciphertexts(ektp, enc_kg, std::get<0>(r1_data.at(i)), randgen);
        }

        // Generate a random unique ID for the group
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(100000, 999999);
        std::string request_id = "request_" + std::to_string(dis(gen));

        std::cout << "Signature generation request created with Request ID: " << request_id << std::endl;
        round1_data.insert(std::make_pair(request_id, std::make_tuple(group_id, enc_kg, kg_dot_G, message_hash_hex_str)));

        return request_id;
    }


    void upload_vn_r1_data(const std::string &request_id, std::vector<std::tuple<BICYCL::CL_HSMqk::CipherText, curve_point>> r1_data, ThreadPool *pool)
    {
        if (is_setup == false)
        {
            throw std::runtime_error("Access Control Contract not setup");
        }

        if (r1_data.size() == 0)
        {
            throw std::runtime_error("No round 1 data to compute on.");
        }

        const ecdsa_curve *curve = get_curve_by_name(SECP256K1_NAME)->params;
        BICYCL::RandGen randgen;
        BICYCL::CL_HSMqk::PublicKey ek = get_public_key(std::get<0>(round1_data.at(request_id)));

        // Determine the segment size and number of threads
        size_t num_threads = std::thread::hardware_concurrency();
        size_t segment_size = (r1_data.size() + num_threads - 1) / num_threads;
        
        std::vector<std::tuple<BICYCL::CL_HSMqk::CipherText, curve_point>> thread_result;
        thread_result.reserve(num_threads);

        // assigning default values
        BICYCL::CL_HSMqk::CipherText default_ciphertext = pp.encrypt(pp.keygen(pp.keygen(randgen)), BICYCL::CL_HSMqk::ClearText(pp, randgen), randgen);
        curve_point default_point;

        for (size_t i = 0; i < num_threads; ++i)
        {
            thread_result.push_back(std::make_tuple(default_ciphertext, default_point));
        }

        // Enqueue tasks to the threadpool
        for (size_t i = 0; i < num_threads; ++i) {
            size_t start_index = i * segment_size;
            if (start_index >= r1_data.size()) break;
            size_t end_index = std::min(start_index + segment_size, r1_data.size());
        
            pool->enqueue(&AccessControlContract::add_r1_data, this, std::ref(r1_data), start_index, end_index, ek, std::ref(thread_result[i])); 
        }

        // Wait for all tasks to complete
        pool->wait();

        // Combine the results from the threads
        curve_point kv_dot_G = std::get<1>(thread_result.at(0));
        BICYCL::CL_HSMqk::CipherText enc_kv = std::get<0>(thread_result.at(0));

        for (size_t i = 1; i < num_threads; ++i) {
            point_add(curve, &std::get<1>(thread_result.at(i)), &kv_dot_G);
            enc_kv = pp.add_ciphertexts(ek, enc_kv, std::get<0>(thread_result.at(i)), randgen);
        }

        std::cout << "ACC then adds this value with group's round 1 values.\n";
        auto &item = round1_data.at(request_id);
        point_add(curve, &kv_dot_G, &std::get<2>(item));
        std::get<1>(round1_data.at(request_id)) = pp.add_ciphertexts(ek, enc_kv, std::get<1>(item), randgen);

        std::cout << "Finally, ACC computes k.G and enc(k) after aggregation.\n";

        // print the final point
        std::cout << "r value (k.G's x coordinate) of the signature: ";
        CryptoUtils::print_bignum256(&std::get<2>(item).x);
    }

    // check for active signature generation requests and return group ids.
    std::vector<std::string> get_active_requests()
    {
        std::vector<std::string> keys;
        for (const auto &item : round1_data)
        {
            keys.push_back(item.first);
        }
        return keys;
    }

    // method to get round1 data given a request id
    std::tuple<std::string, BICYCL::CL_HSMqk::CipherText, curve_point, std::string> get_round1_data(const std::string &request_id)
    {
        if (round1_data.find(request_id) != round1_data.end())
        {
            return round1_data.at(request_id);
        }
        else
        {
            throw std::runtime_error("Request ID not found"); // Throw an error if the request ID is not found
        }
    }

    void add_r2_data(std::vector<std::tuple<BICYCL::CL_HSMqk::CipherText, BICYCL::CL_HSMqk::CipherText, BICYCL::CL_HSMqk::CipherText>> &r2_data, 
                     size_t start_index, 
                     size_t end_index, 
                     BICYCL::CL_HSMqk::PublicKey ektp,
                     std::tuple<BICYCL::CL_HSMqk::CipherText, BICYCL::CL_HSMqk::CipherText, BICYCL::CL_HSMqk::CipherText> &result) {
        BICYCL::RandGen randgen;

        BICYCL::CL_HSMqk::CipherText enc_p = std::get<0>(r2_data.at(start_index));
        BICYCL::CL_HSMqk::CipherText enc_p_times_x = std::get<1>(r2_data.at(start_index));
        BICYCL::CL_HSMqk::CipherText enc_p_times_k = std::get<2>(r2_data.at(start_index));

        for (size_t i = start_index + 1; i < end_index; ++i)
        {
            enc_p = pp.add_ciphertexts(ektp, enc_p, std::get<0>(r2_data.at(i)), randgen);
            enc_p_times_x = pp.add_ciphertexts(ektp, enc_p_times_x, std::get<1>(r2_data.at(i)), randgen);
            enc_p_times_k = pp.add_ciphertexts(ektp, enc_p_times_k, std::get<2>(r2_data.at(i)), randgen);
        }

        result = std::make_tuple(enc_p, enc_p_times_x, enc_p_times_k);
    }

    void upload_r2_data(std::string request_id, std::vector<std::tuple<BICYCL::CL_HSMqk::CipherText, BICYCL::CL_HSMqk::CipherText, BICYCL::CL_HSMqk::CipherText>> &r2_data, ThreadPool *pool)
    {
        if (is_setup == false)
        {
            throw std::runtime_error("Access Control Contract not setup");
        }
        BICYCL::RandGen randgen;

        std::cout << "ACC combines round 2 data of validator network and the user group.\n";
        std::string group_id = std::get<0>(get_round1_data(request_id));
        BICYCL::CL_HSMqk::PublicKey ektp = get_public_key(group_id);

        // Determine the segment size and number of threads
        size_t num_threads = std::thread::hardware_concurrency();
        size_t segment_size = (r2_data.size() + num_threads - 1) / num_threads;
        
        std::vector<std::tuple<BICYCL::CL_HSMqk::CipherText, BICYCL::CL_HSMqk::CipherText, BICYCL::CL_HSMqk::CipherText>> thread_result;
        thread_result.reserve(num_threads);

        // assigning default values
        BICYCL::CL_HSMqk::CipherText default_ciphertext = pp.encrypt(pp.keygen(pp.keygen(randgen)), BICYCL::CL_HSMqk::ClearText(pp, randgen), randgen);

        for (size_t i = 0; i < num_threads; ++i)
        {
            thread_result.push_back(std::make_tuple(default_ciphertext, default_ciphertext, default_ciphertext));
        }

        // Enqueue tasks to the threadpool
        for (size_t i = 0; i < num_threads; ++i) {
            size_t start_index = i * segment_size;
            if (start_index >= r2_data.size()) break;
            size_t end_index = std::min(start_index + segment_size, r2_data.size());
        
            pool->enqueue(&AccessControlContract::add_r2_data, this, std::ref(r2_data), start_index, end_index, ektp, std::ref(thread_result[i])); 
        }

        // Wait for all tasks to complete
        pool->wait();

        // Combine the results from the threads
        BICYCL::CL_HSMqk::CipherText enc_p = std::get<0>(thread_result.at(0));
        BICYCL::CL_HSMqk::CipherText enc_p_times_x = std::get<1>(thread_result.at(0));
        BICYCL::CL_HSMqk::CipherText enc_p_times_k = std::get<2>(thread_result.at(0));

        for (size_t i = 1; i < num_threads; ++i) {
            enc_p = pp.add_ciphertexts(ektp, enc_p, std::get<0>(thread_result[i]), randgen);
            enc_p_times_x = pp.add_ciphertexts(ektp, enc_p_times_x, std::get<1>(thread_result[i]), randgen);
            enc_p_times_k = pp.add_ciphertexts(ektp, enc_p_times_k, std::get<2>(thread_result[i]), randgen);
        }

        round2_data.insert(std::make_pair(request_id, std::make_tuple(group_id, enc_p, enc_p_times_x, enc_p_times_k)));
    }

    std::tuple<std::string, BICYCL::CL_HSMqk::CipherText, BICYCL::CL_HSMqk::CipherText, BICYCL::CL_HSMqk::CipherText> get_round2_data(const std::string &request_id)
    {
        if (round2_data.find(request_id) != round2_data.end())
        {
            return round2_data.at(request_id);
        }
        else
        {
            throw std::runtime_error("Request ID not found"); // Throw an error if the request ID is not found
        }
    }

    void nucomp_parallel(const BICYCL::QFI decryptions[],
                         size_t start_index,
                         size_t end_index,
                         BICYCL::QFI &result)
    {
        result = decryptions[start_index];
        for (size_t i = start_index + 1; i < end_index; ++i)
        {
            pp.Cl_Delta().nucomp(result, result, decryptions[i]);
        }
    }

    BICYCL::CL_HSMqk::ClearText agg_partial_decryptions_parallel(const BICYCL::QFI decryptions[], 
                                                                 const size_t n, 
                                                                 const BICYCL::CL_HSMqk::CipherText &c,
                                                                 ThreadPool *pool)
    {
        BICYCL::QFI c2 = c.c2();

        // Determine the segment size and number of threads
        size_t num_threads = std::thread::hardware_concurrency();
        size_t segment_size = (n + num_threads - 1) / num_threads;

        BICYCL::QFI thread_result[num_threads];

        // Enqueue tasks to the threadpool
        for (size_t i = 0; i < num_threads; ++i) {
            size_t start_index = i * segment_size;
            if (start_index >= n) break;
            size_t end_index = std::min(start_index + segment_size, n);
        
            pool->enqueue(&AccessControlContract::nucomp_parallel, this, decryptions, start_index, end_index, std::ref(thread_result[i]));
        }

        // Wait for all tasks to complete
        pool->wait();
        BICYCL::QFI temp = thread_result[0];

        for (size_t i = 1; i < num_threads; ++i)
        {
            pp.Cl_Delta().nucomp(temp, temp, thread_result[i]);
        }

        pp.Cl_Delta().nucompinv(c2, c2, temp);
        return BICYCL::CL_HSMqk::ClearText (pp, pp.dlog_in_F(c2));
    }

    void upload_r3_data(std::string request_id, BICYCL::QFI partial_decryptions[], size_t size, ThreadPool *pool)
    {
        if (is_setup == false)
        {
            throw std::runtime_error("Access Control Contract not setup");
        }

        std::cout << "ACC combines round 3 data of validator network and the user group.\n";

        std::tuple<std::string, BICYCL::CL_HSMqk::CipherText, curve_point, std::string> r1_data = get_round1_data(request_id);
        std::tuple<std::string, BICYCL::CL_HSMqk::CipherText, BICYCL::CL_HSMqk::CipherText, BICYCL::CL_HSMqk::CipherText> r2_data = get_round2_data(request_id);

        BICYCL::CL_HSMqk::CipherText enc_p_times_k = std::get<3>(r2_data);
        BICYCL::CL_HSMqk::ClearText w = agg_partial_decryptions_parallel(partial_decryptions, size, enc_p_times_k, pool);

        std::cout << "ACC computes encrypted signature from the encrypted data.\n";
        std::string group_id = std::get<0>(r1_data);
        BICYCL::CL_HSMqk::PublicKey ektp = get_public_key(group_id);
        BICYCL::RandGen randgen;

        // contract then computes w^-1
        BICYCL::Mpz w_inv(0UL);
        BICYCL::Mpz::mod_inverse(w_inv, w, pp.q());

        // message hash
        BICYCL::Mpz message_hash(0UL);
        std::string message_hash_str = std::get<3>(r1_data);
        BICYCL::Mpz::set_str(message_hash, message_hash_str.c_str(), 16);

        // contract computes w^-1 * hash(m) and w^-1 * r (used in computing encrypted signature)
        BICYCL::Mpz w_inv_times_message_hash(0UL);
        BICYCL::Mpz::mul(w_inv_times_message_hash, w_inv, message_hash);
        // BICYCL::Mpz::mod(w_inv_times_message_hash, w_inv_times_message_hash, q);

        // r (x-coordinate of k.G)
        BICYCL::Mpz r(0UL);
        uint8_t r_bytes[32];
        bn_write_be(&(std::get<2>(r1_data)).x, r_bytes);
        BICYCL::Mpz::from_bytes(r, r_bytes, 32);

        BICYCL::Mpz w_inv_times_r(0UL);
        BICYCL::Mpz::mul(w_inv_times_r, w_inv, r);
        // BICYCL::Mpz::mod(w_inv_times_r, w_inv_times_r, q);

        BICYCL::CL_HSMqk::CipherText enc_s = pp.add_ciphertexts(ektp,
                                                                pp.scal_ciphertexts(ektp, std::get<1>(r2_data), w_inv_times_message_hash, randgen),
                                                                pp.scal_ciphertexts(ektp, std::get<2>(r2_data), w_inv_times_r, randgen),
                                                                randgen);

        round3_data.insert(std::make_pair(request_id, std::make_tuple(group_id, enc_s)));
    }

    std::tuple<std::string, BICYCL::CL_HSMqk::CipherText> get_round3_data(const std::string &request_id)
    {
        if (round3_data.find(request_id) != round3_data.end())
        {
            return round3_data.at(request_id);
        }
        else
        {
            throw std::runtime_error("Request ID not found"); // Throw an error if the request ID is not found
        }
    }

    void upload_r4_data(std::string request_id, BICYCL::QFI partial_decryptions[], size_t size)
    {
        if (is_setup == false)
        {
            throw std::runtime_error("Access Control Contract not setup");
        }

        std::cout << "ACC combines round 4 data of validator network and the user group.\n";

        std::tuple<std::string, BICYCL::CL_HSMqk::CipherText> r3_data = get_round3_data(request_id);

        BICYCL::CL_HSMqk::CipherText enc_s = std::get<1>(r3_data);
        BICYCL::CL_HSMqk::ClearText s = pp.agg_partial_decryptions(partial_decryptions, size, enc_s);

        std::string group_id = std::get<0>(r3_data);

        std::cout << "ACC stores the final decrypted signature.\n";
        round4_data.insert(std::make_pair(request_id, std::make_tuple(group_id, s)));
    }

    std::tuple<std::string, BICYCL::CL_HSMqk::ClearText> get_round4_data(const std::string &request_id)
    {
        if (round4_data.find(request_id) != round4_data.end())
        {
            return round4_data.at(request_id);
        }
        else
        {
            throw std::runtime_error("Request ID not found"); // Throw an error if the request ID is not found
        }
    }
};

class Validator
{
private:
    BICYCL::Mpz decryption_key_share; // Private data member
    bool is_setup = false;

public:
    // Setter for decryption_key_share
    void set_decryption_key_share(BICYCL::Mpz keyShare)
    {
        decryption_key_share = keyShare;
        is_setup = true;
    }

    void compute_r1_data(BICYCL::CL_HSMqk pp,
                         const BICYCL::CL_HSMqk::PublicKey ek,
                         std::tuple<BICYCL::CL_HSMqk::CipherText, curve_point> &result_data,
                         std::mutex *m)
    {
        BICYCL::RandGen randgen;

        const ecdsa_curve *curve = get_curve_by_name(SECP256K1_NAME)->params;
        curve_point kvi_dot_G;

        bignum256 kvi_bn;
        CryptoUtils::rand_bignum256(&curve->order, &kvi_bn);

        uint8_t kvi_bytes[32];
        bn_write_be(&kvi_bn, kvi_bytes);

        BICYCL::Mpz kvi(0UL);
        BICYCL::Mpz::from_bytes(kvi, kvi_bytes, 32);

        BICYCL::CL_HSMqk::ClearText kvi_cleartext = BICYCL::CL_HSMqk::ClearText(pp, kvi);
        BICYCL::CL_HSMqk::CipherText enc_kvi = pp.encrypt(ek, kvi_cleartext, randgen);

        std::lock_guard<std::mutex> lock(*m);
        scalar_multiply(curve, &kvi_bn, &kvi_dot_G);

        result_data = std::make_tuple(enc_kvi, kvi_dot_G);
    }

    void compute_r2_data(BICYCL::CL_HSMqk pp,
                         const BICYCL::CL_HSMqk::PublicKey ek,
                         BICYCL::CL_HSMqk::CipherText enc_k,
                         BICYCL::CL_HSMqk::CipherText enc_x,
                         std::tuple<BICYCL::CL_HSMqk::CipherText, BICYCL::CL_HSMqk::CipherText, BICYCL::CL_HSMqk::CipherText> &result_data)
    {
        BICYCL::RandGen randgen;
        BICYCL::Mpz pvi = randgen.random_mpz(pp.q());

        BICYCL::CL_HSMqk::ClearText pvi_cleartext = BICYCL::CL_HSMqk::ClearText(pp, pvi);
        BICYCL::CL_HSMqk::CipherText enc_pvi = pp.encrypt(ek, pvi_cleartext, randgen);

        BICYCL::CL_HSMqk::CipherText enc_pvi_times_x = pp.scal_ciphertexts(ek, enc_x, pvi, randgen);
        BICYCL::CL_HSMqk::CipherText enc_pvi_times_k = pp.scal_ciphertexts(ek, enc_k, pvi, randgen);

        result_data = std::make_tuple(enc_pvi, enc_pvi_times_x, enc_pvi_times_k);
    }

    void compute_r3_data(BICYCL::CL_HSMqk pp,
                         BICYCL::CL_HSMqk::CipherText enc_p_times_k,
                         BICYCL::QFI &part_dec)
    {
        BICYCL::Mpz lambda(0UL);
        lambda = string("1");

        part_dec = pp.partial_decrypt(BICYCL::CL_HSMqk::SecretKey(pp, decryption_key_share), lambda, enc_p_times_k);
    }

    void compute_r4_data(BICYCL::CL_HSMqk pp,
                         BICYCL::CL_HSMqk::CipherText enc_s,
                         BICYCL::QFI &part_dec)
    {
        BICYCL::Mpz lambda(0UL);
        lambda = string("1");

        part_dec = pp.partial_decrypt(BICYCL::CL_HSMqk::SecretKey(pp, decryption_key_share), lambda, enc_s);
    }
};

class ValidatorNetwork
{
private:
    std::vector<Validator> validators; // Vector to store Validator objects
    AccessControlContract *acc;
    BICYCL::CL_HSMqk pp; // CL_HSMqk object to store public parameters
    bool is_setup = false;

public:
    // Constructor that initializes the vector with a specific number of Validators
    ValidatorNetwork(int num_validators, AccessControlContract *acc, BICYCL::CL_HSMqk pp) : validators(num_validators), acc(acc), pp(pp)
    {
        std::cout << "Setting up validator network with " << num_validators << " validators\n";
    }

    size_t size()
    {
        return validators.size();
    }

    // method to mock validator network's DKG, assigns decryption key shares to validators
    // and returns the network's public encryption key 'ekv'. Currently, all validators are
    // required to reconstruct the secret decryption key.
    void mock_dkg_and_setup_acc()
    {
        if (is_setup)
        {
            throw std::runtime_error("Validator network already setup");
        }

        std::cout << "Performing DKG in the validator network...\n";

        BICYCL::RandGen randgen;

        // adding 1 to account for user group's decryption keys.
        int k = std::ceil(std::log2(validators.size())) + 1;

        BICYCL::Mpz decryption_key_share_bound(0UL);
        BICYCL::Mpz::divby2k(decryption_key_share_bound, pp.secretkey_bound(), k);

        BICYCL::Mpz skv_raw(0UL);
        for (Validator &validator : validators)
        {
            BICYCL::Mpz skv_i = BICYCL::Mpz(randgen.random_mpz(decryption_key_share_bound));
            BICYCL::Mpz::add(skv_raw, skv_raw, skv_i);
            validator.set_decryption_key_share(skv_i);
        }

        BICYCL::CL_HSMqk::SecretKey skv = BICYCL::CL_HSMqk::SecretKey(pp, skv_raw);
        BICYCL::CL_HSMqk::PublicKey ekv = pp.keygen(skv);

        acc->set_encryption_key(ekv);
        is_setup = true;
    }

    std::vector<std::string> check_sig_requests()
    {
        return acc->get_active_requests();
    }

    void participate_in_round1(std::string request_id, ThreadPool *pool)
    {
        std::cout << "Validator network retrieves round 1 data and public encryption key from ACC.\n";
        // get round 1 data
        std::tuple<std::string, BICYCL::CL_HSMqk::CipherText, curve_point, std::string> round1_data = acc->get_round1_data(request_id);
        // get encryption key from group id
        BICYCL::CL_HSMqk::PublicKey ek = acc->get_public_key(std::get<0>(round1_data));

        std::vector<std::tuple<BICYCL::CL_HSMqk::CipherText, curve_point>> validator_network_r1_data;

        BICYCL::RandGen randgen;
        BICYCL::CL_HSMqk::CipherText default_ciphertext = pp.encrypt(pp.keygen(pp.keygen(randgen)), BICYCL::CL_HSMqk::ClearText(pp, randgen), randgen);
        curve_point default_point;

        for (size_t i = 0; i < validators.size(); ++i)
        {
            validator_network_r1_data.push_back(std::make_tuple(default_ciphertext, default_point));
        }

        std::cout << "Validators in the validator network compute round 1 data parallelly and upload it to ACC.\n";

        std::mutex data_mutex;
        size_t index = 0;

        for (Validator &validator : validators)
        {
            pool->enqueue(&Validator::compute_r1_data, &validator, pp, ek, std::ref(validator_network_r1_data[index]), &data_mutex);
            index += 1;
        }

        pool->wait();

        return acc->upload_vn_r1_data(request_id, validator_network_r1_data, pool);
    }

    void participate_in_round2(std::string request_id, std::vector<std::tuple<BICYCL::CL_HSMqk::CipherText, BICYCL::CL_HSMqk::CipherText, BICYCL::CL_HSMqk::CipherText>> &round2_data, size_t offset, ThreadPool *pool)
    {
        std::cout << "Validator Network retrieves round 1 data, encrypted signing key and public encryption key from ACC.\n";
        // get round 1 data
        std::tuple<std::string, BICYCL::CL_HSMqk::CipherText, curve_point, std::string> round1_data = acc->get_round1_data(request_id);
        // get encrypted signing key from group id
        BICYCL::CL_HSMqk::CipherText enc_x = acc->get_cipher_text(std::get<0>(round1_data));
        // get encryption key from group id
        BICYCL::CL_HSMqk::PublicKey ek = acc->get_public_key(std::get<0>(round1_data));

        BICYCL::CL_HSMqk::CipherText enc_k = std::get<1>(round1_data);

        std::cout << "Validators in the validator network compute round 2 data parallelly.\n";
        size_t index = offset;
        for (Validator &validator : validators)
        {
            pool->enqueue(&Validator::compute_r2_data, &validator, pp, ek, enc_k, enc_x, std::ref(round2_data[index]));
            ++index;
        }
    }

    void participate_in_round3(std::string request_id, BICYCL::QFI partial_decryptions[], size_t offset, ThreadPool *pool)
    {
        std::cout << "Validator Network retrieves round 2 data from ACC.\n";
        // get round 2 data
        std::tuple<std::string, BICYCL::CL_HSMqk::CipherText, BICYCL::CL_HSMqk::CipherText, BICYCL::CL_HSMqk::CipherText> round2_data = acc->get_round2_data(request_id);

        std::cout << "Validators in the validator network compute round 3 data parallelly.\n";
        size_t index = offset;
        for (Validator &validator : validators)
        {
            pool->enqueue(&Validator::compute_r3_data, &validator, pp, std::get<3>(round2_data), std::ref(partial_decryptions[index]));
            ++index;
        }
    }

    void participate_in_round4(std::string request_id, BICYCL::QFI partial_decryptions[], size_t offset, ThreadPool *pool)
    {
        std::cout << "Validator Network retrieves round 3 data from ACC.\n";
        // get round 3 data
        std::tuple<std::string, BICYCL::CL_HSMqk::CipherText> round3_data = acc->get_round3_data(request_id);

        std::cout << "Validators in the validator network compute round 4 data parallelly.\n";
        size_t index = offset;
        for (Validator &validator : validators)
        {
            pool->enqueue(&Validator::compute_r4_data, &validator, pp, std::get<1>(round3_data), std::ref(partial_decryptions[index]));
            ++index;
        }
    }
};

class User
{
private:
    uint32_t id;                      // starts from 1
    BICYCL::Mpz decryption_key_share; // Private data member
    bool is_setup = false;

public:
    // Setter for decryption_key_share
    void set_decryption_key_share(uint32_t uid, BICYCL::Mpz keyShare)
    {
        id = uid;
        decryption_key_share = keyShare;
        is_setup = true;
    }

    void compute_r1_data(BICYCL::CL_HSMqk &pp,
                         const BICYCL::CL_HSMqk::PublicKey &ek,
                         std::vector<std::tuple<BICYCL::CL_HSMqk::CipherText, curve_point>> &result_data,
                         std::mutex *m)
    {
        BICYCL::RandGen randgen;
        curve_point kgi_dot_G;
        const ecdsa_curve *curve = get_curve_by_name(SECP256K1_NAME)->params;

        bignum256 kgi_bn;
        CryptoUtils::rand_bignum256(&curve->order, &kgi_bn);

        uint8_t kgi_bytes[32];
        bn_write_be(&kgi_bn, kgi_bytes);

        BICYCL::Mpz kgi(0UL);
        BICYCL::Mpz::from_bytes(kgi, kgi_bytes, 32);

        BICYCL::CL_HSMqk::ClearText kgi_cleartext = BICYCL::CL_HSMqk::ClearText(pp, kgi);
        BICYCL::CL_HSMqk::CipherText enc_kgi = pp.encrypt(ek, kgi_cleartext, randgen);

        std::lock_guard<std::mutex> lock(*m);
        scalar_multiply(curve, &kgi_bn, &kgi_dot_G);

        result_data.push_back(std::make_tuple(enc_kgi, kgi_dot_G)); 
    }


    void compute_r2_data(BICYCL::CL_HSMqk pp,
                         const BICYCL::CL_HSMqk::PublicKey ek,
                         BICYCL::CL_HSMqk::CipherText enc_k,
                         BICYCL::CL_HSMqk::CipherText enc_x,
                         std::tuple<BICYCL::CL_HSMqk::CipherText, BICYCL::CL_HSMqk::CipherText, BICYCL::CL_HSMqk::CipherText> &result_data)
    {

        BICYCL::RandGen randgen;
        BICYCL::Mpz pgi = randgen.random_mpz(pp.q());

        BICYCL::CL_HSMqk::ClearText pgi_cleartext = BICYCL::CL_HSMqk::ClearText(pp, pgi);
        BICYCL::CL_HSMqk::CipherText enc_pgi = pp.encrypt(ek, pgi_cleartext, randgen);

        BICYCL::CL_HSMqk::CipherText enc_pgi_times_x = pp.scal_ciphertexts(ek, enc_x, pgi, randgen);
        BICYCL::CL_HSMqk::CipherText enc_pgi_times_k = pp.scal_ciphertexts(ek, enc_k, pgi, randgen);

        result_data = std::make_tuple(enc_pgi, enc_pgi_times_x, enc_pgi_times_k);
    }

    void compute_r3_data(BICYCL::CL_HSMqk pp,
                         BICYCL::CL_HSMqk::CipherText enc_p_times_k,
                         BICYCL::QFI &part_dec)
    {
        BICYCL::Mpz lambda(0UL);
        lambda = string("1");

        part_dec = pp.partial_decrypt(BICYCL::CL_HSMqk::SecretKey(pp, decryption_key_share), lambda, enc_p_times_k);
    }

    void compute_r4_data(BICYCL::CL_HSMqk pp,
                         BICYCL::CL_HSMqk::CipherText enc_s,
                         BICYCL::QFI &part_dec)
    {
        BICYCL::Mpz lambda(0UL);
        lambda = string("1");

        part_dec = pp.partial_decrypt(BICYCL::CL_HSMqk::SecretKey(pp, decryption_key_share), lambda, enc_s);
    }
};

class UserGroup
{
private:
    std::vector<User> users; // Vector to store User objects
    PolicyContract *policy;
    AccessControlContract *acc;
    BICYCL::CL_HSMqk::PublicKey encryption_key;
    BICYCL::CL_HSMqk pp; // CL_HSMqk object to store public parameters
    bool is_setup = false;
    bool is_registered = false;
    std::string group_id;

public:
    // Constructor that initializes the vector with a specific number of Users
    UserGroup(int num_users, PolicyContract *policy, AccessControlContract *acc, BICYCL::CL_HSMqk pp) : users(num_users), policy(policy), acc(acc),
                                                                                                        encryption_key(BICYCL::CL_HSMqk::PublicKey(
                                                                                                            pp, BICYCL::CL_HSMqk::SecretKey(pp, (BICYCL::RandGen()).random_mpz(pp.secretkey_bound())))),
                                                                                                        pp(pp)
    {
        std::cout << "Setting up User Group with " << num_users << " group members.\n";
    }

    size_t size()
    {
        return users.size();
    }

    // method to mock user group's DKG, assigns decryption key shares to users
    // and returns the group's public encryption key 'ekg'. Currently, all users are
    // required to reconstruct the secret decryption key. Note that in future,
    // multiple keys will be distributed because of multiple thresholds.
    std::string mock_dkg_and_register_group()
    {
        if (is_setup)
        {
            throw std::runtime_error("User group already setup");
        }

        BICYCL::RandGen randgen;

        // adding 1 to account for validator network's decryption keys.
        int k = std::ceil(std::log2(users.size())) + 1;

        BICYCL::Mpz decryption_key_share_bound(0UL);
        BICYCL::Mpz::divby2k(decryption_key_share_bound, pp.secretkey_bound(), k);

        BICYCL::Mpz skg_raw(0UL);

        std::cout << "Perform User Group DKG for decryption key...\n";
        uint32_t i = 1;
        for (User &user : users)
        {
            BICYCL::Mpz skg_i = BICYCL::Mpz(randgen.random_mpz(decryption_key_share_bound));
            BICYCL::Mpz::add(skg_raw, skg_raw, skg_i);
            user.set_decryption_key_share(i, skg_i);
            i += 1;
        }

        BICYCL::CL_HSMqk::SecretKey skg = BICYCL::CL_HSMqk::SecretKey(pp, skg_raw);
        BICYCL::CL_HSMqk::PublicKey ekg = pp.keygen(skg);

        BICYCL::CL_HSMqk::PublicKey ekv = acc->get_encryption_key();

        std::cout << "Combining validator network's encryption key with group's encryption key.\n";

        BICYCL::QFI ekg_raw = ekg.elt();
        BICYCL::QFI ekv_raw = ekv.elt();

        BICYCL::QFI ek_raw;

        // TODO - Check if Cl_Delta is okay for all cases (compact variant)
        pp.Cl_Delta().nucomp(ek_raw, ekg_raw, ekv_raw);

        BICYCL::CL_HSMqk::PublicKey ek = BICYCL::CL_HSMqk::PublicKey(pp, ek_raw);
        encryption_key = ek;

        // Mocking the generation of random group public key and encrypted signing key
        const ecdsa_curve *curve = get_curve_by_name(SECP256K1_NAME)->params;

        std::cout << "Generating encrypted signing key for the group...\n";

        // private signing key (never actually exists in plain form)
        bignum256 x_bn;
        CryptoUtils::rand_bignum256(&curve->order, &x_bn);
        uint8_t signing_private_key_bytes[32];
        bn_write_be(&x_bn, signing_private_key_bytes);

        BICYCL::Mpz x(0UL);
        BICYCL::Mpz::from_bytes(x, signing_private_key_bytes, 32);

        uint8_t group_pub_key[65];
        ecdsa_get_public_key65(curve, signing_private_key_bytes, group_pub_key);
        
        // print group_pub_key
        std::cout << "Group verification key: ";
        CryptoUtils::print_bytes(group_pub_key, 65);

        // encrypted signing key is publicly available
        BICYCL::CL_HSMqk::ClearText x_cleartext = BICYCL::CL_HSMqk::ClearText(pp, x);
        BICYCL::CL_HSMqk::CipherText enc_x = pp.encrypt(ek, x_cleartext, randgen); // encrypted signing key

        group_id = acc->register_user_group(policy, ek, enc_x, group_pub_key);
        is_setup = true;

        return group_id;
    }

    std::string initiate_signature_generation(std::string message_hash_hex_str, ThreadPool *pool)
    {
        if (!is_setup)
        {
            throw std::runtime_error("User group not setup");
        }

        // print message
        std::cout << "Message to sign: " << message_hash_hex_str << std::endl;

        // for this poc
        uint32_t threshold = users.size();

        // hash of the message to be signed.
        BICYCL::Mpz message_hash(0UL);
        BICYCL::Mpz::set_str(message_hash, message_hash_hex_str.c_str(), 16);

        uint8_t message_hash_bytes[32];
        BICYCL::Mpz::to_bytes(message_hash, message_hash_bytes, 32); // it causes a problem if starting bytes are 00

        std::vector<std::tuple<BICYCL::CL_HSMqk::CipherText, curve_point>> user_group_r1_data;
        user_group_r1_data.reserve(users.size());

        std::cout << "Users in the user group compute round 1 data parallelly before sending signature generation request to ACC.\n";

        std::mutex data_mutex;
        size_t index = 0;
        for (User &user : users)
        {
            pool->enqueue(&User::compute_r1_data, &user, pp, encryption_key, std::ref(user_group_r1_data), &data_mutex);
            ++index;
        }

        pool->wait();

        std::cout << "Uploading the request parameters (t and M) along with round 1 data to ACC.\n";
        return acc->upload_ug_r1_data(group_id, threshold, message_hash_hex_str, user_group_r1_data);
    }

    void participate_in_round2(std::string request_id, std::vector<std::tuple<BICYCL::CL_HSMqk::CipherText, BICYCL::CL_HSMqk::CipherText, BICYCL::CL_HSMqk::CipherText>> &round2_data, size_t offset, ThreadPool *pool)
    {
        // get round 1 data
        std::cout << "User group retrieves round 1 data and encrypted signing key from ACC.\n";
        std::tuple<std::string, BICYCL::CL_HSMqk::CipherText, curve_point, std::string> round1_data = acc->get_round1_data(request_id);

        // get encrypted signing key
        BICYCL::CL_HSMqk::CipherText enc_x = acc->get_cipher_text(std::get<0>(round1_data));
        BICYCL::CL_HSMqk::CipherText enc_k = std::get<1>(round1_data);

        std::cout << "Users in the user group compute round 2 data parallelly.\n";
        size_t index = offset;
        for (User &user : users)
        {
            pool->enqueue(&User::compute_r2_data, &user, pp, encryption_key, enc_k, enc_x, std::ref(round2_data[index]));
            ++index;
        }
    }

    void participate_in_round3(std::string request_id, BICYCL::QFI partial_decryptions[], size_t offset, ThreadPool *pool)
    {
        // get round 2 data
        std::cout << "User group retrieves round 2 data from ACC.\n";
        std::tuple<std::string, BICYCL::CL_HSMqk::CipherText, BICYCL::CL_HSMqk::CipherText, BICYCL::CL_HSMqk::CipherText> round2_data = acc->get_round2_data(request_id);

        std::cout << "Users in the user group compute round 3 data parallelly.\n";
        size_t index = offset;
        for (User &user : users)
        {
            pool->enqueue(&User::compute_r3_data, &user, pp, std::get<3>(round2_data), std::ref(partial_decryptions[index]));
            ++index;
        }
    }

    void participate_in_round4(std::string request_id, BICYCL::QFI partial_decryptions[], size_t offset, ThreadPool *pool)
    {
        // get round 3 data
        std::cout << "User group retrieves round 3 data from ACC.\n";
        std::tuple<std::string, BICYCL::CL_HSMqk::CipherText> round3_data = acc->get_round3_data(request_id);

        std::cout << "Users in the user group compute round 4 data parallelly.\n";
        size_t index = offset;
        for (User &user : users)
        {
            pool->enqueue(&User::compute_r4_data, &user, pp, std::get<1>(round3_data), std::ref(partial_decryptions[index]));
            ++index;
        }
    }
};

int main(int argc, char* argv[])
{
    // Check if the correct number of arguments are passed
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <number_of_users> <number_of_validators> <security_level>\n";
        return 1;
    }

    // Convert command line arguments from strings to integers
    int number_of_users = std::atoi(argv[1]);
    int number_of_validators = std::atoi(argv[2]);
    int security_level = std::atoi(argv[3]);

    // Validate security level
    if (security_level != 128 && security_level != 256) {
        std::cerr << "Error: Security level must be either 128 or 256.\n";
        return 1;
    }

    size_t k = 1;
    const ecdsa_curve *curve = get_curve_by_name(SECP256K1_NAME)->params;
    BICYCL::RandGen randgen;

    ThreadPool pool(std::thread::hardware_concurrency());
    std::cout << "Created ThreadPool of size: " << std::thread::hardware_concurrency() << "\n";

    // Setting up the the public parameters where message space is the curve order of SECP256K1
    BICYCL::Mpz q(0UL);
    q = string("115792089237316195423570985008687907852837564279074904382605163141518161494337");

    BICYCL::Mpz p(0UL);

    if (security_level == 128) {
        p = string("43193454325736827482201750544140481404008634846532639450547372755026427405991826700723249528255049592325421795981143653914129884011306525120950744025035685855512695631736186909872842154602977956593529422235466196411674948283663780511726036070197880528511138509668014828860898402959171280613034578545771859496996085810790266951981368946435378708249190772155224039947314675049652489452429836293254671938998647344790956467563810567741334084883797963767053599118688017207673143");
    } else if (security_level == 256) {
        p = string ("18162659420476506244413899955873056039189126922056845480032142867443044062547737609595400815586415693442107398687192520033132150181492794933559538297169597362513634417261752318910589620230655148591000545739393647489358969239828175446943072123909611576032188421712348810334719982487055378242959377856456215580786926934463584071028076071697683772543753807565840321970401447542241241440740152712738757651049133498513361341544091878946024159470118454826395344604747651882311618739796455773166441059612341833624528960159357526099782595726203287558033173472552072322172418326412899878599880998801656108307944262446105942907727187448886677694724879883196438685160111217043051893930552040177062272443085934772801004322558940176706738218072343625205182388345074274711495024428597770830733258073816860067136150841876903078873109532159417650394957993730726916008083158428712673522649614680398837928784169489054950302445228225796652290286004105294592749781604616876668706031722924705454246527806652373891211331111561836142813415753871453254546219902033827618368259391980900615511617687273265288771889191101025935249962339676750387432829515461932986082221170815113578523503210282985388722188002200489238868628364715568641769286094590701848803043593464313512088224481009878875601271864226519775365382549278086391433047919878981001384713552785810526234467984926199913445309334969179294983198805783660234577808758212261681158552057515664328931933487376751755258239249690772173778395073782451548897885708535209703585975795907058251800456387679810899906500271428284578707852013386078542402930654702400845624245878806838375059891857419154842582118569686130930273522084293521704081743508162833028585321947555706174712379336870496920297557603");
    }

    std::cout << "Setting up class group...\n";
    BICYCL::CL_HSMqk pp = BICYCL::CL_HSMqk(q, k, p);
    std::cout << "Class group successfully setup.\n";

    AccessControlContract acc = AccessControlContract(pp);

    std::cout << "\n";
    std::cout << "===========================================================================================================================================================\n";
    std::cout << "Setting up the validator network - Distributed Key Generation\n";
    std::cout << "===========================================================================================================================================================\n";
    auto start = std::chrono::high_resolution_clock::now();

    ValidatorNetwork validator_network = ValidatorNetwork(number_of_validators, &acc, pp);
    validator_network.mock_dkg_and_setup_acc();

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    std::cout << "Validator network setup complete.\n";
    std::cout << "===========================================================================================================================================================\n";
    std::cout << "Time taken to setup validator network: " << duration.count() << " milliseconds" << std::endl;
    std::cout << "===========================================================================================================================================================\n";

    std::cout << "\n";
    std::cout << "===========================================================================================================================================================\n";
    std::cout << "Setting up the user group - Distributed Key Generation and group registration\n";
    std::cout << "===========================================================================================================================================================\n";
    start = std::chrono::high_resolution_clock::now();

    GroupsPolicyContract groups_policy = GroupsPolicyContract();
    UserGroup user_group(number_of_users, &groups_policy, &acc, pp);

    std::string group_id = user_group.mock_dkg_and_register_group();

    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    std::cout << "User group setup and registration complete.\n";
    std::cout << "===========================================================================================================================================================\n";
    std::cout << "Time taken to setup and register user group: " << duration.count() << " milliseconds" << std::endl;
    std::cout << "===========================================================================================================================================================\n";

    std::cout << "\n";
    std::cout << "===========================================================================================================================================================\n";
    std::cout << "Signature Generation Process\n";
    std::cout << "===========================================================================================================================================================\n";
    auto sig_gen_start = std::chrono::high_resolution_clock::now();

    std::cout << "Round - 1\n";
    std::cout << "===========================================================================================================================================================\n";
    start = std::chrono::high_resolution_clock::now();

    std::string message_hash_hex_str = "01d4043e3e023fef92f6546a06a4cfbb60b938ec145d4809706a138061fea544";
    std::string ug_req_id = user_group.initiate_signature_generation(message_hash_hex_str, &pool);

    std::cout << "Validator network checks for active signature generation requests.\n";
    std::vector<std::string> request_ids = validator_network.check_sig_requests();

    // print all group ids
    std::cout << "Found the following active signature generation requests: \n";
    for (std::string request_id : request_ids)
    {
        std::cout << request_id << std::endl;
    }

    std::cout << "Validator network participates in signature generation for the request: " << request_ids.at(0) << std::endl;
    validator_network.participate_in_round1(request_ids.at(0), &pool);

    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    std::cout << "===========================================================================================================================================================\n";
    std::cout << "Time taken in Round - 1: " << duration.count() << " milliseconds" << std::endl;
    std::cout << "===========================================================================================================================================================\n";

    std::cout << "Round - 2\n";
    std::cout << "===========================================================================================================================================================\n";
    start = std::chrono::high_resolution_clock::now();

    // Round 2 data to broadcast - enc(pi), enc(pi * x), enc(pi * k)
    std::cout << "Setting up storage for round 2 data...\n";
    std::vector<std::tuple<BICYCL::CL_HSMqk::CipherText, BICYCL::CL_HSMqk::CipherText, BICYCL::CL_HSMqk::CipherText>> round2_data;

    // assigning default values
    BICYCL::CL_HSMqk::CipherText default_ciphertext = pp.encrypt(pp.keygen(pp.keygen(randgen)), BICYCL::CL_HSMqk::ClearText(pp, randgen), randgen);

    for (size_t i = 0; i < user_group.size() + validator_network.size(); ++i)
    {
        round2_data.push_back(std::make_tuple(default_ciphertext, default_ciphertext, default_ciphertext));
    }

    std::cout << "User group participates in round 2.\n";
    user_group.participate_in_round2(ug_req_id, round2_data, 0, &pool);

    std::cout << "Validator network participates in round 2 parallelly.\n";
    validator_network.participate_in_round2(request_ids.at(0), round2_data, user_group.size(), &pool);

    // wait for all threads to finish
    pool.wait();

    std::cout << "Round 2 data is broadcasted to ACC.\n";
    acc.upload_r2_data(ug_req_id, round2_data, &pool);

    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    std::cout << "===========================================================================================================================================================\n";
    std::cout << "Time taken in Round - 2: " << duration.count() << " milliseconds" << std::endl;
    std::cout << "===========================================================================================================================================================\n";

    std::cout << "Round - 3\n";
    std::cout << "===========================================================================================================================================================\n";
    start = std::chrono::high_resolution_clock::now();

    // Round 3 data to broadcast - part_dec(enc(p * k))
    std::cout << "Setting up storage for round 3 data...\n";
    BICYCL::QFI w_partial_decryptions[user_group.size() + validator_network.size()];

    std::cout << "User group participates in round 3.\n";
    user_group.participate_in_round3(ug_req_id, w_partial_decryptions, 0, &pool);

    std::cout << "Validator network participates in round 3 parallelly.\n";
    validator_network.participate_in_round3(request_ids.at(0), w_partial_decryptions, user_group.size(), &pool);

    // wait for all threads to finish
    pool.wait();

    std::cout << "Round 3 data is broadcasted to ACC.\n";
    acc.upload_r3_data(ug_req_id, w_partial_decryptions, user_group.size() + validator_network.size(), &pool);

    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    std::cout << "===========================================================================================================================================================\n";
    std::cout << "Time taken in Round - 3: " << duration.count() << " milliseconds" << std::endl;
    std::cout << "===========================================================================================================================================================\n";

    std::cout << "Round - 4\n";
    std::cout << "===========================================================================================================================================================\n";
    start = std::chrono::high_resolution_clock::now();

    // Round 4 data to broadcast - part_dec(enc_s)
    std::cout << "Setting up storage for round 4 data...\n";
    BICYCL::QFI s_partial_decryptions[user_group.size() + validator_network.size()];

    std::cout << "User group participates in round 4.\n";
    user_group.participate_in_round4(ug_req_id, s_partial_decryptions, 0, &pool);

    std::cout << "Validator network participates in round 4 parallelly.\n";
    validator_network.participate_in_round4(request_ids.at(0), s_partial_decryptions, user_group.size(), &pool);

    // wait for all threads to finish
    pool.wait();

    std::cout << "Round 4 data is broadcasted to ACC.\n";
    acc.upload_r4_data(ug_req_id, s_partial_decryptions, user_group.size() + validator_network.size());

    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    std::cout << "Signature Generated." << std::endl;

    std::cout << "===========================================================================================================================================================\n";
    std::cout << "Time taken in Round - 4: " << duration.count() << " milliseconds" << std::endl;
    std::cout << "===========================================================================================================================================================\n";

    auto sig_gen_stop = std::chrono::high_resolution_clock::now();
    auto sig_gen_duration = std::chrono::duration_cast<std::chrono::milliseconds>(sig_gen_stop - sig_gen_start);
    std::cout << "Total time taken to generate signature: " << sig_gen_duration.count() << " milliseconds" << std::endl;
    std::cout << "===========================================================================================================================================================\n";

    auto r1_data = acc.get_round1_data(ug_req_id);
    auto r4_data = acc.get_round4_data(ug_req_id);

    uint8_t sig_bytes[64];

    // copy the x-coordinate of k.G to the signature
    curve_point k_dot_G = std::get<2>(r1_data);
    bn_write_be(&k_dot_G.x, sig_bytes);
    BICYCL::Mpz::to_bytes(std::get<1>(r4_data), sig_bytes + 32, 32);

    // print signature
    std::cout << "\n";
    std::cout << "Signature: ";
    CryptoUtils::print_bytes(sig_bytes, 64);

    uint8_t verification_key_bytes[65];
    CryptoUtils::hex_string_to_byte_array(acc.get_verification_key(group_id).c_str(), verification_key_bytes, 65);

    uint8_t message_hash_bytes[32];
    CryptoUtils::hex_string_to_byte_array(std::get<3>(r1_data).c_str(), message_hash_bytes, 32);

    if (ecdsa_verify_digest(curve, verification_key_bytes, sig_bytes, message_hash_bytes) == 0)
    {
        printf("Signature verified successfully\n");
    }
    else
    {
        printf("Signature verification failed\n");
    }

    return 0;
}
