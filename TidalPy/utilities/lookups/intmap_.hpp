#pragma once

#include <vector>
#include <algorithm>
#include <cstddef>
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

        if (this->data.empty() || key_ref > this->data.back().first.reference)
        {
            this->data.emplace_back(key, value);
            return;
        }

        const std::size_t index = this->p_lower_index(key_ref);
        if (index < this->data.size() && this->data[index].first.reference == key_ref)
        {
            this->data[index].second = value;
        } else
        {
            this->data.insert(this->data.begin() + index, {key, value});
        }
    }

    ValueType get(bool& o_found, const KeyType& key) const
    {
        o_found = true;

        const std::size_t index = this->p_find_index(key.reference);
        if (index < this->data.size())
        {
            return this->data[index].second;
        }

        o_found = false;
        return ValueType(); // Default construct (0.0 or 0j)
    }

    const ValueType* get_ptr(bool& o_found, const KeyType& key) const
    {
        o_found = true;

        const std::size_t index = this->p_find_index(key.reference);
        if (index < this->data.size())
        {
            return &(this->data[index].second);  // Points into the vector; a later set() can invalidate it.
        }

        o_found = false;
        return nullptr;
    }

    auto begin() { return this->data.begin(); }
    auto end() { return this->data.end(); }

    auto begin() const { return this->data.begin(); }
    auto end() const { return this->data.end(); }

private:
    // Index of the first entry whose key is not below key_ref; data is kept sorted by the packed key reference.
    std::size_t p_lower_index(RefKeyType key_ref) const
    {
        const auto it = std::lower_bound(this->data.begin(), this->data.end(), key_ref,
            [](const std::pair<KeyType, ValueType>& entry, RefKeyType k) { return entry.first.reference < k; });
        return static_cast<std::size_t>(it - this->data.begin());
    }

    // Index of the entry holding key_ref, or data.size() when there is none.
    std::size_t p_find_index(RefKeyType key_ref) const
    {
        const std::size_t index = this->p_lower_index(key_ref);
        if (index < this->data.size() && this->data[index].first.reference == key_ref)
        {
            return index;
        }
        return this->data.size();
    }
};
