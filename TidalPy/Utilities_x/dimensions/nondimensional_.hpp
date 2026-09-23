#pragma once
/* Unit-conversion scales for non-dimensionalized solves, following Martens (2016, PhD Thesis, CalTech,
 * DOI: 10.7907/Z9N29TX7, ~p. 99): the time scale from 1/(pi G rho_bulk) (frequency independent), the
 * length scale from the mean radius, the density scale from the bulk density; mass and pascal follow.
 */

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
        double mean_radius,
        double bulk_density
    )
    {
        this->second2_conversion = 1. / (TidalPyConstants::d_PI * c_get_G() * bulk_density);
        this->second_conversion  = std::sqrt(this->second2_conversion);
        this->length_conversion  = mean_radius;
        this->length3_conversion = mean_radius * mean_radius * mean_radius;
        this->density_conversion = bulk_density;
        this->mass_conversion    = bulk_density * this->length3_conversion;
        this->pascal_conversion  = \
            this->mass_conversion / (this->length_conversion * this->second2_conversion);
    }
};
