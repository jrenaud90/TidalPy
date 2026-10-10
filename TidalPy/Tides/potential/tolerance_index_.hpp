#pragma once
/*
 * tolerance_index_.hpp - c_ToleranceIndex: the earliest stored value a relative-tolerance match accepts, found without
 * scanning every stored value.
 *
 * Entries carry a tag, which must match exactly, and a value, which must pass the caller's match predicate. A lookup
 * returns what a scan in insertion order returns: the insertion index of the first entry that matches. Finite values
 * sit in a list sorted by (tag, value), and a finite query checks only the stretch of it within
 * [v (1 - rtol), v / (1 - rtol)], widened a little beyond the rounding of either end. That window holds every value a
 * relative tolerance of rtol can accept, so the predicate, not the window, decides. Entries the window cannot place (a
 * non-finite value) are checked on every lookup, a non-finite query checks every entry, and an rtol outside [0, 0.5]
 * turns the window off, each in insertion order through the predicate as well.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

template <typename Predicate>
class c_ToleranceIndex {
public:
    // predicate(stored, query) decides a match; rtol is the relative tolerance it applies.
    c_ToleranceIndex(double rtol, Predicate predicate) :
        p_rtol(rtol),
        p_predicate(predicate),
        p_windowed((rtol >= 0.0) && (rtol <= 0.5)) {}

    std::size_t size() const noexcept { return this->p_entries.size(); }

    // The insertion index of the first entry with this tag whose value matches, or -1 when none does.
    std::ptrdiff_t find(std::int64_t tag, double value) const {
        std::ptrdiff_t first = -1;
        this->for_each_match(tag, value, [&](std::size_t index) {
            if ((first < 0) || (static_cast<std::ptrdiff_t>(index) < first)) {
                first = static_cast<std::ptrdiff_t>(index);
            }
        });
        return first;
    }

    // visit(index) for every entry with this tag whose value matches, in no particular order.
    template <typename Visit>
    void for_each_match(std::int64_t tag, double value, const Visit& visit) const {
        if (!this->p_windowed || !std::isfinite(value)) {
            for (std::size_t index = 0; index < this->p_entries.size(); ++index) {
                const c_Entry& entry = this->p_entries[index];
                if ((entry.tag == tag) && this->p_predicate(entry.value, value)) { visit(index); }
            }
            return;
        }
        // Past any rounding of the window ends, relative for normal numbers and absolute near zero.
        const double margin = std::abs(value) * 1.0e-12 + 16.0 * std::numeric_limits<double>::denorm_min();
        const double near_end = value * (1.0 - this->p_rtol);
        const double far_end = value / (1.0 - this->p_rtol);
        const double low = std::min(near_end, far_end) - margin;
        const double high = std::max(near_end, far_end) + margin;
        auto it = std::lower_bound(
            this->p_sorted.begin(), this->p_sorted.end(), c_Entry{tag, low, 0}, c_EntryOrder());
        for (; it != this->p_sorted.end(); ++it) {
            if ((it->tag != tag) || (it->value > high)) { break; }
            if (this->p_predicate(it->value, value)) { visit(it->index); }
        }
        for (const c_Entry& entry : this->p_unsorted) {
            if ((entry.tag == tag) && this->p_predicate(entry.value, value)) { visit(entry.index); }
        }
    }

    // Adds an entry and returns its insertion index.
    std::size_t insert(std::int64_t tag, double value) {
        const c_Entry entry{tag, value, this->p_entries.size()};
        this->p_entries.push_back(entry);
        if (!std::isfinite(value)) {
            this->p_unsorted.push_back(entry);
        } else if (this->p_windowed) {
            this->p_sorted.insert(
                std::upper_bound(this->p_sorted.begin(), this->p_sorted.end(), entry, c_EntryOrder()), entry);
        }
        return entry.index;
    }

private:
    struct c_Entry {
        std::int64_t tag;
        double value;
        std::size_t index;
    };

    struct c_EntryOrder {
        bool operator()(const c_Entry& a, const c_Entry& b) const noexcept {
            return (a.tag < b.tag) || ((a.tag == b.tag) && (a.value < b.value));
        }
    };

    double p_rtol;
    Predicate p_predicate;
    bool p_windowed;
    std::vector<c_Entry> p_entries;    // every entry, in insertion order
    std::vector<c_Entry> p_sorted;     // the finite entries by (tag, value), when windowed
    std::vector<c_Entry> p_unsorted;   // the non-finite entries
};
