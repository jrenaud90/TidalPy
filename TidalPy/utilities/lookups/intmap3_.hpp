#pragma once

#include <vector>
#include <algorithm>
#include <cstdint>
#include <utility>
#include <complex> // Required for complex support

// Define KeyType
typedef uint64_t KeyType3;

// Correct OFFSET (int32_t) to prevent overflow
const int32_t OFFSET = 32768;

// Helper to pack key
inline KeyType3 convet_3key(int16_t l, int16_t m, int16_t p) {
    return (static_cast<KeyType3>(l + OFFSET) << 32) | 
           (static_cast<KeyType3>(m + OFFSET) << 16) | 
           static_cast<KeyType3>(p + OFFSET);
}


template <typename T>
class c_IntMap3 {
private:
    std::vector<std::pair<KeyType3, T>> data;

public:
    c_IntMap3() {
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

    // Changed 'double' to 'T'
    void set(int16_t l, int16_t m, int16_t p, T value) {
        KeyType3 key = convet_3key(l, m, p);

        if (data.empty() || key > data.back().first) {
            data.emplace_back(key, value);
            return;
        }

        auto it = std::lower_bound(data.begin(), data.end(), key, 
            [](const auto& entry, KeyType3 k) { return entry.first < k; });

        if (it != data.end() && it->first == key) {
            it->second = value;
        } else {
            data.insert(it, {key, value});
        }
    }

    // Changed 'double' to 'T'
    // Note: We return T by value here for simplicity, or we can use reference param
    T get(bool& o_found, int16_t l, int16_t m, int16_t p) const {
        KeyType3 key = convet_3key(l, m, p);
        o_found = true;
        
        auto it = std::lower_bound(data.begin(), data.end(), key, 
            [](const auto& entry, KeyType3 k) { return entry.first < k; });
            
        if (it != data.end() && it->first == key) {
            return it->second;
        }

        o_found = false;
        return T(); // Default construct (0.0 or 0j)
    }
};