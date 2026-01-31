#include "intmap4_.hpp"

// CONSTANTS FOR PACKING
// We support signed 16-bit integers: -32,768 to +32,767
const int32_t OFFSET = 32768;

inline KeyType4 convet_4key(int16_t l, int16_t m, int16_t p, int16_t q)
{
    // Helper to pack the key
    return (static_cast<KeyType4>(l + OFFSET) << 48) | 
           (static_cast<KeyType4>(m + OFFSET) << 32) | 
           (static_cast<KeyType4>(p + OFFSET) << 16) | 
           static_cast<KeyType4>(q + OFFSET);
}

c_IntMap4::c_IntMap4()
{
    // Reserve capacity
    data.reserve(35);
}

void c_IntMap4::reserve(size_t n)
{
    // Reserve memory if you know the rough size (prevents reallocation)
    this->data.reserve(n);
}

void c_IntMap4::set(int16_t l, int16_t m, int16_t p, int16_t q, double value)
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
}

size_t c_IntMap4::size() const
{
    return this->data.size();
}

double c_IntMap4::get(bool& o_found, int16_t l, int16_t m, int16_t p, int16_t q) const
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
}

const std::vector<std::pair<KeyType4, double>>& c_IntMap4::get_raw_data() const
{
    // Since we are usually iterating iterate in order (2,0,0,0), (2,0,1,0)...
    // Expose the raw vector for iteration. 
    // This is faster than calling get() in a loop.
    return this->data;
}

void c_IntMap4::clear()
{
    this->data.clear();
}
