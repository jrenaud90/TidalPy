#pragma once

#include <cmath>
#include <utility>

#include "intmap_.hpp"
#include "keys_.hpp"

typedef std::pair<c_IntMap<c_Key3, double>, c_IntMap<c_Key2, c_IntMap<c_Key1, double>>> EccentricityFuncOutput;
