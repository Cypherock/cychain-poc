#pragma once

#include <string>

class PolicyContract {
public:
    virtual ~PolicyContract() {}

    // Pure virtual function that must be implemented by any concrete class deriving from this interface
    virtual bool check_compliance(uint32_t threshold, std::string message_hash_hex_str) = 0;
};
