#include "intmap3_.hpp"

inline KeyType convet_3key(size_t l, size_t m, size_t p)
{
    // Helper to pack the key
    return (static_cast<KeyType>(l) << 20) | 
           (static_cast<KeyType>(m) << 10) | 
           static_cast<KeyType>(p);
}

c_IntMap3::c_IntMap3()
{
    // Reserve capacity
    data.reserve(35);
}

void c_IntMap3::reserve(size_t n)
{
    // Reserve memory if you know the rough size (prevents reallocation)
    this->data.reserve(n);
}

void c_IntMap3::set(size_t l, size_t m, size_t p, double value)
{
    // Insert or Update
    // O(N) worst case, but O(1) if you insert in order (append).
    KeyType key = convet_3key(l, m, p);

    // Optimization: If inserting in strict order, just push_back
    if (this->data.empty() || key > this->data.back().first) {
        this->data.emplace_back(key, value);
        return;
    }

    // Otherwise, find position to keep it sorted
    auto it = std::lower_bound(this->data.begin(), this->data.end(), key, 
        [](const auto& entry, KeyType k) { return entry.first < k; });

    if (it != this->data.end() && it->first == key) {
        it->second = value; // Update existing
    } else {
        this->data.insert(it, {key, value}); // Insert new
    }
}

size_t c_IntMap3::size() const
{
    return this->data.size();
}

double c_IntMap3::get(bool& o_found, size_t l, size_t m, size_t p) const
{
    // Lookup
    // O(log N) - Binary Search. For N=200, this is ~8 comparisons
    KeyType key = convet_3key(l, m, p);
    o_found = true;
    
    auto it = std::lower_bound(this->data.begin(), this->data.end(), key, 
        [](const auto& entry, KeyType k) { return entry.first < k; });
        
    if (it != this->data.end() && it->first == key) {
        return it->second;
    }

    // Not found. Generally we want this to return 0.0 in this case.
    o_found = false;
    return 0.0;
}

const std::vector<std::pair<KeyType, double>>& c_IntMap3::get_raw_data() const
{
    // Since we are usually iterating iterate in order (2,0,0), (2,0,1)...
    // Expose the raw vector for iteration. 
    // This is faster than calling get() in a loop.
    return this->data;
}

void c_IntMap3::clear()
{
    this->data.clear();
}
