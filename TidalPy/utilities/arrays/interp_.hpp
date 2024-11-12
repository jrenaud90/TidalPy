#pragma once

#include <limits>

static const double EPS = std::numeric_limits<double>::epsilon();

size_t cf_binary_search_with_guess(double key, double* array, size_t length, size_t guess);
