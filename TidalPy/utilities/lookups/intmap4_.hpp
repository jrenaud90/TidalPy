#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdint>

typedef uint64_t KeyType4;

// Helper to pack the key
inline KeyType4 convet_4key(int16_t l, int16_t m, int16_t p, int16_t q);

class c_IntMap4 {
private:
    // We pack (l, m, p) into a single 32-bit integer.
    // Assuming l, m, p < 65,535. We allocate 16 bits per number (up to 65,535).
    // Layout: [16 bits Unused | 16 bits m | 16 bits p | 16 bits q]
    
    // Store data contiguously.
    // pair.first = Packed Key, pair.second = Value
    std::vector<std::pair<KeyType4, double>> data;
public:
    c_IntMap4();

    // Reserve memory if you know the rough size (prevents reallocation)
    void reserve(size_t n);

    // Insert or Update
    // O(N) worst case, but O(1) if you insert in order (append).
    void set(int16_t l, int16_t m, int16_t p, int16_t q, double value);

    // Lookup
    // O(log N) - Binary Search. For N=200, this is ~8 comparisons.
    double get(bool& o_found, int16_t l, int16_t m, int16_t p, int16_t q) const;
    size_t size() const;

    // Since we are usually iterating iterate in order (2,0,0), (2,0,1)...
    // Expose the raw vector for iteration. 
    // This is faster than calling get() in a loop.
    const std::vector<std::pair<KeyType4, double>>& get_raw_data() const;
    
    // Clear the vector
    void clear();
};