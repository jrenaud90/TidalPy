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

    c_IntMap()
    {
        this->data.reserve(30);
    }
    
    c_IntMap(size_t n)
    {
        this->data.reserve(n);
    }

    void reserve(size_t n)
    {
        this->data.reserve(n);
    }

    void clear()
    {
        this->data.clear();
    }
    
    size_t size() const
    {
        return this->data.size();
    }

    void set(const KeyType& key, const ValueType& value)
    {

        // The 64-bit packed signature of this key.
        RefKeyType key_ref = key.reference;

        if (data.empty() || key_ref > data.back().first.reference)
        {
            data.emplace_back(key, value);
            return;
        }

        auto it = std::lower_bound(data.begin(), data.end(), key_ref,
            [](const auto& entry, RefKeyType k) { return entry.first.reference < k; });


        if (it != data.end() && it->first.reference == key_ref)
        {
            it->second = value;
        } else
        {
            data.insert(it, {key, value});
        }
    }

    ValueType get(bool& o_found, const KeyType& key) const
    {
        o_found = true;

        RefKeyType key_ref = key.reference;
        
        auto it = std::lower_bound(data.begin(), data.end(), key_ref, 
            [](const auto& entry, RefKeyType k) { return entry.first.reference < k; });
            
        if (it != data.end() && it->first.reference == key_ref)
        {
            return it->second;
        }

        o_found = false;
        return ValueType(); // Default construct (0.0 or 0j)
    }

    const ValueType* get_ptr(bool& o_found, const KeyType& key) const
    {
        o_found = true;

        RefKeyType key_ref = key.reference;
        
        auto it = std::lower_bound(data.begin(), data.end(), key_ref, 
            [](const auto& entry, RefKeyType k) { return entry.first.reference < k; });
            
        if (it != data.end() && it->first.reference == key_ref)
        {
            return &(it->second);  // Points into the vector; a later set() can invalidate it.
        }

        return nullptr;
    }

    auto begin() { return data.begin(); }
    auto end() { return data.end(); }

    auto begin() const { return data.begin(); }
    auto end() const { return data.end(); }
};
