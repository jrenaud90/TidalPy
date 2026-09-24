#pragma once
/* Model-name handling shared by the physics-model factories, which match names and aliases case-insensitively. */

#include <algorithm>
#include <cctype>
#include <string>

namespace tidalpy {

// ASCII lower case, the form every model name and alias is compared in.
inline std::string c_to_lower(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(),
                   [](unsigned char character) { return static_cast<char>(std::tolower(character)); });
    return text;
}

} // namespace tidalpy
