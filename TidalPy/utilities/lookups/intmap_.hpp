#pragma once

#include <vector>
#include <algorithm>
#include <cstdint>
#include <utility>
#include <tuple>
#include <complex>

#include "keys_.hpp"


template <typename KeyType, typename ValueType>
class c_IntMap {

public:
    std::vector<std::pair<KeyType, ValueType>> data;

    c_IntMap() {
        data.reserve(35);
    }

    void reserve(size_t n) {
        data.reserve(n);
    }

    void clear() {
        data.clear();
    }
    
    size_t size() const {
        return data.size();
    }

    void set(KeyType& key, ValueType value) {

        // Get 64-bit integer that stores a signature of this key.
        RefKeyType key_ref = key.reference;

        if (data.empty() || key_ref > data.back().first.reference) {
            data.emplace_back(key, value);
            return;
        }

        auto it = std::lower_bound(data.begin(), data.end(), key_ref,
            [](const auto& entry, RefKeyType k) { return entry.first.reference < k; }); // lambda function that finds the first element that does not compare less than key


        if (it != data.end() && it->first.reference == key_ref) {
            // Update existing entry.
            it->second = value;
        } else {
            data.insert(it, {key, value});
        }
    }

    // Note: We return T by value here for simplicity, or we can use reference param
    ValueType get(bool& o_found, KeyType& key) const {
        o_found = true;

        // Get 64-bit integer that stores a signature of this key.
        RefKeyType key_ref = key.reference;
        
        auto it = std::lower_bound(data.begin(), data.end(), key_ref, 
            [](const auto& entry, RefKeyType k) { return entry.first.reference < k; });
            
        if (it != data.end() && it->first.reference == key_ref) {
            return it->second;
        }

        o_found = false;
        return ValueType(); // Default construct (0.0 or 0j)
    }
};
