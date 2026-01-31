#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdint>

typedef uint32_t KeyType;

// Helper to pack the key
inline KeyType convet_3key(size_t l, size_t m, size_t p);

class c_IntMap3 {
private:
    // We pack (l, m, p) into a single 32-bit integer.
    // Assuming l, m, p < 100. We allocate 10 bits per number (up to 1023).
    // Layout: [ 2 bits unused | 10 bits l | 10 bits m | 10 bits p ]
    
    // Store data contiguously.
    // pair.first = Packed Key, pair.second = Value
    std::vector<std::pair<KeyType, double>> data;
public:
    c_IntMap3();

    // Reserve memory if you know the rough size (prevents reallocation)
    void reserve(size_t n);

    // Insert or Update
    // O(N) worst case, but O(1) if you insert in order (append).
    void set(size_t l, size_t m, size_t p, double value);

    // Lookup
    // O(log N) - Binary Search. For N=200, this is ~8 comparisons.
    double get(bool& o_found, size_t l, size_t m, size_t p) const;
    size_t size() const;

    // Since we are usually iterating iterate in order (2,0,0), (2,0,1)...
    // Expose the raw vector for iteration. 
    // This is faster than calling get() in a loop.
    const std::vector<std::pair<KeyType, double>>& get_raw_data() const;
    
    // Clear the vector
    void clear();
};