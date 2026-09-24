#pragma once
/* Unit-conversion scales for non-dimensionalized solves, following Martens (2016, PhD Thesis, CalTech,
 * DOI: 10.7907/Z9N29TX7, ~p. 99): the time scale from 1/(pi G rho_bulk) (frequency independent), the
 * length scale from the mean radius, the density scale from the bulk density; mass and pascal follow.
 */

#include <cmath>
#include "../../constants_.hpp"

// Every scale is NaN until built from a planet's radius and density. Copyable and assignable.
class c_NonDimensionalScales
{
public:
    double second2_conversion = TidalPyConstants::d_NAN;
    double second_conversion  = TidalPyConstants::d_NAN;
    double length_conversion  = TidalPyConstants::d_NAN;
    double length3_conversion = TidalPyConstants::d_NAN;
    double density_conversion = TidalPyConstants::d_NAN;
    double mass_conversion    = TidalPyConstants::d_NAN;
    double pascal_conversion  = TidalPyConstants::d_NAN;

    c_NonDimensionalScales() = default;
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
