#pragma once

#include <cmath>
#include "../../constants_.hpp"

class c_NonDimensionalScales
{
public:
    double second2_conversion;
    double second_conversion;
    double length_conversion;
    double length3_conversion;
    double density_conversion;
    double mass_conversion;
    double pascal_conversion;

    c_NonDimensionalScales() = default;
    ~c_NonDimensionalScales() = default;
    c_NonDimensionalScales(const c_NonDimensionalScales&) = default;
    c_NonDimensionalScales(c_NonDimensionalScales&&) = default;
    c_NonDimensionalScales(
        double frequency,
        double mean_radius,
        double bulk_density
    )
    {
        this->second2_conversion = 1. / (TidalPyConstants::d_PI * tidalpy_config_ptr->d_G * bulk_density);
        this->second_conversion  = std::sqrt(this->second2_conversion);
        this->length_conversion  = mean_radius;
        this->length3_conversion = mean_radius * mean_radius * mean_radius;
        this->density_conversion = bulk_density;
        this->mass_conversion    = bulk_density * this->length3_conversion;
        this->pascal_conversion  = \
            this->mass_conversion / (this->length_conversion * this->second2_conversion);
    }
};
