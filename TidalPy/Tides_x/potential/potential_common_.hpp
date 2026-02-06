#pragma once

#include <cstdint>
#include <vector>

#include "intmap_.hpp"
#include "keys_.hpp"
#include "numerics_.hpp"
#include "constants_.hpp"

struct FrequencyStorage
{
    size_t num_instances;
    double frequency;

    FrequencyStorage(double freq) :
        num_instances(1),
        frequency(freq)
    {
    }
};

typedef c_IntMap<c_Key4, double> ModeMap;
typedef c_IntMap<c_Key4, size_t> UniqueFreqIndexMap;
typedef std::vector<FrequencyStorage> UniqueFreqMap;


void record_unique_frequencies(
        c_Key4& lmpq_key,
        double frequency,
        UniqueFreqIndexMap& frequency_index_map,
        UniqueFreqMap& frequency_map)
{
    // TODO: Do we want to keep zero frequencies?
    if (c_isclose(frequency, 0.0, 1.0e-9, tidalpy_config_ptr->d_MIN_FREQUENCY))
    {
        /* code */
    }
    
    // Checks if this frequency has been recorded yet. If not then it adds it to the frequency map. 
    bool found = false;
    FrequencyStorage* frequency_storage_ptr = nullptr;
    for (size_t i = 0; i < frequency_map.size(); i++)
    {
        frequency_storage_ptr = &frequency_map[i];
        if (c_isclose(frequency, frequency_storage_ptr->frequency, 1.0e-9, 0))
        {   
            // Increment the number of times this frequency has shown up.
            frequency_storage_ptr->num_instances += 1;
            // Record that this mode corresponds to this unique frequency.
            frequency_index_map.set(lmpq_key, i);
            found = true;
            // Break out of the frequency loop.
            break;
        }
    }

    if (!found)
    {
        // Did not find the frequency, so this is a unique one!
        frequency_index_map.set(lmpq_key, frequency_map.size());
        frequency_map.emplace_back(frequency);
    }
}
