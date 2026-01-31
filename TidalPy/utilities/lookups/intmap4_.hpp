#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdint>

typedef uint64_t KeyType4;

// Helper to pack the key
// We support signed 16-bit integers: -32,768 to +32,767
const int32_t OFFSET = 32768;

inline KeyType4 convet_4key(int16_t l, int16_t m, int16_t p, int16_t q)
{
    // Helper to pack the key
    return (static_cast<KeyType4>(l + OFFSET) << 48) | 
           (static_cast<KeyType4>(m + OFFSET) << 32) | 
           (static_cast<KeyType4>(p + OFFSET) << 16) | 
           static_cast<KeyType4>(q + OFFSET);
};

template <typename T>
class c_IntMap4 {
private:
    // We pack (l, m, p) into a single 32-bit integer.
    // Assuming l, m, p < 65,535. We allocate 16 bits per number (up to 65,535).
    // Layout: [16 bits Unused | 16 bits m | 16 bits p | 16 bits q]
    
    // Store data contiguously.
    // pair.first = Packed Key, pair.second = Value
    std::vector<std::pair<KeyType4, T>> data;
public:
    c_IntMap4()
    {
        // Reserve capacity
        data.reserve(35);
    };

    void reserve(size_t n)
    {
        // Reserve memory if you know the rough size (prevents reallocation)
        this->data.reserve(n);
    };

    void set(int16_t l, int16_t m, int16_t p, int16_t q, T value)
    {
        // Insert or Update
        // O(N) worst case, but O(1) if you insert in order (append).
        KeyType4 key = convet_4key(l, m, p, q);

        // Optimization: If inserting in strict order, just push_back
        if (this->data.empty() || key > this->data.back().first) {
            this->data.emplace_back(key, value);
            return;
        }

        // Otherwise, find position to keep it sorted
        auto it = std::lower_bound(this->data.begin(), this->data.end(), key, 
            [](const auto& entry, KeyType4 k) { return entry.first < k; });

        if (it != this->data.end() && it->first == key) {
            it->second = value; // Update existing
        } else {
            this->data.insert(it, {key, value}); // Insert new
        }
    };

    T get(bool& o_found, int16_t l, int16_t m, int16_t p, int16_t q) const
    {
        // Lookup
        // O(log N) - Binary Search. For N=200, this is ~8 comparisons
        KeyType4 key = convet_4key(l, m, p, q);
        o_found = true;
        
        auto it = std::lower_bound(this->data.begin(), this->data.end(), key, 
            [](const auto& entry, KeyType4 k) { return entry.first < k; });
            
        if (it != this->data.end() && it->first == key) {
            return it->second;
        }

        // Not found. Generally we want this to return 0.0 in this case.
        o_found = false;
        return 0.0;
    };

    size_t size() const
    {
        return this->data.size();
    };

    const std::vector<std::pair<KeyType4, T>>& get_raw_data() const
    {
        // Since we are usually iterating iterate in order (2,0,0,0), (2,0,1,0)...
        // Expose the raw vector for iteration. 
        // This is faster than calling get() in a loop.
        return this->data;
    };
    
    // Clear the vector
    void clear()
    {
        this->data.clear();
    };
};