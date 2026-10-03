#pragma once
/* Name handling shared by the physics-model factories and the parameter specs: case-insensitive matching, and the
 * closest accepted name for an error message. */

#include <algorithm>
#include <cctype>
#include <string>
#include <vector>

namespace tidalpy {

// ASCII lower case, the form every model name and alias is compared in.
inline std::string c_to_lower(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(),
                   [](unsigned char character) { return static_cast<char>(std::tolower(character)); });
    return text;
}

// Levenshtein distance between two strings, case-insensitive.
inline std::size_t c_edit_distance(const std::string& first, const std::string& second) {
    const std::string a = c_to_lower(first);
    const std::string b = c_to_lower(second);
    std::vector<std::size_t> previous(b.size() + 1);
    std::vector<std::size_t> current(b.size() + 1);
    for (std::size_t j = 0; j <= b.size(); ++j) { previous[j] = j; }
    for (std::size_t i = 1; i <= a.size(); ++i) {
        current[0] = i;
        for (std::size_t j = 1; j <= b.size(); ++j) {
            const std::size_t substitution = previous[j - 1] + ((a[i - 1] == b[j - 1]) ? 0 : 1);
            current[j] = std::min({previous[j] + 1, current[j - 1] + 1, substitution});
        }
        std::swap(previous, current);
    }
    return previous[b.size()];
}

// The candidate closest to `name`, or an empty string when none is close enough to be a likely misspelling: within a
// third of the longer name's length, or one that `name` begins (a key missing its unit suffix, say).
inline std::string c_closest_name(const std::string& name, const std::vector<std::string>& candidates) {
    const std::string lowered = c_to_lower(name);
    std::string best;
    std::size_t best_distance = static_cast<std::size_t>(-1);
    for (const std::string& candidate : candidates) {
        const std::string candidate_lowered = c_to_lower(candidate);
        const bool prefix = !lowered.empty() && (candidate_lowered.rfind(lowered, 0) == 0);
        const std::size_t distance = prefix ? 0 : c_edit_distance(lowered, candidate_lowered);
        const std::size_t limit = std::max(lowered.size(), candidate_lowered.size()) / 3;
        if ((distance <= limit) && (distance < best_distance)) {
            best = candidate;
            best_distance = distance;
        }
    }
    return best;
}

// " (did you mean 'x'?)" for the closest candidate, or nothing.
inline std::string c_did_you_mean(const std::string& name, const std::vector<std::string>& candidates) {
    const std::string closest = c_closest_name(name, candidates);
    return closest.empty() ? std::string() : " (did you mean '" + closest + "'?)";
}

} // namespace tidalpy
