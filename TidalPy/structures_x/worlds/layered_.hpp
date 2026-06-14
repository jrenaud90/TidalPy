#pragma once
/*
 * layered_.hpp — c_LayeredWorld: a world built from an ordered stack of layers.
 *
 * Inherits c_BaseWorld. Owns its layers as std::unique_ptr<c_BaseLayer> (inner to
 * outer, index 0 = innermost) and provides whole-planet aggregates (total mass,
 * internal radiogenic heating) and geometry validation. The whole-planet EOS and
 * radial solves (which walk all layers) are added as methods on this class in
 * later phases.
 *
 * Binary format (20-byte header + payload):
 *   header: class_id = BinaryClassID::LayeredWorld (201)
 *   payload: [all c_BaseWorld fields] + layer_count (uint64_t)
 *   Each layer's own complete binary record (with its recursively-serialized
 *   physics sub-models) follows, in index order, as separate appended records.
 */

#include <cstdint>
#include <istream>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <vector>

#include "base_.hpp"
#include "../layers/factory_.hpp"

namespace tidalpy {

// -------------------------------------------------------------------------------
// c_LayeredWorld
// -------------------------------------------------------------------------------
class c_LayeredWorld : public c_BaseWorld {
public:
    // Relative tolerance used when checking layer-boundary continuity.
    static constexpr double d_LAYER_CONTINUITY_RTOL = 1.0e-6;

    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------
    c_LayeredWorld() = default;

    explicit c_LayeredWorld(const c_WorldConfig& cfg) : c_BaseWorld(cfg) {}

    ~c_LayeredWorld() override = default;

    // -----------------------------------------------------------------------
    // Layer ownership
    // -----------------------------------------------------------------------
    // Add a layer, inner to outer. Validates that its inner radius matches the
    // current outermost radius (0 for the first layer) within the relative
    // tolerance. Throws std::invalid_argument on a gap/overlap.
    void add_layer(std::unique_ptr<c_BaseLayer> layer) {
        if (!layer) {
            throw std::invalid_argument("TidalPy: cannot add a null layer to a world");
        }
        const double prev_outer = this->p_layers.empty() ? 0.0
                                : this->p_layers.back()->get_radius_outer();
        const double inner      = layer->get_radius_inner();
        const double tol        = d_LAYER_CONTINUITY_RTOL
                                * (prev_outer > 1.0 ? prev_outer : 1.0);
        if (std::abs(inner - prev_outer) > tol) {
            throw std::invalid_argument(
                "TidalPy: layer geometry is not continuous — inner radius does not "
                "match the previous layer's outer radius (add layers inner-to-outer)");
        }
        this->p_layers.push_back(std::move(layer));
    }

    // True if `layer` would continue the stack: its inner radius matches the
    // current outermost radius (0 for the first layer) within tolerance. Lets a
    // caller check continuity before transferring ownership via add_layer.
    bool accepts_layer(const c_BaseLayer& layer) const noexcept {
        const double prev_outer = this->p_layers.empty() ? 0.0
                                : this->p_layers.back()->get_radius_outer();
        const double tol        = d_LAYER_CONTINUITY_RTOL
                                * (prev_outer > 1.0 ? prev_outer : 1.0);
        return std::abs(layer.get_radius_inner() - prev_outer) <= tol;
    }

    // Non-owning observer pointer to the layer at index (throws out_of_range).
    c_BaseLayer* get_layer(std::size_t index) const {
        if (index >= this->p_layers.size()) {
            throw std::out_of_range("TidalPy: layer index out of range");
        }
        return this->p_layers[index].get();
    }

    std::size_t get_num_layers() const noexcept { return this->p_layers.size(); }

    // -----------------------------------------------------------------------
    // Whole-planet aggregates (const, MKS)
    // -----------------------------------------------------------------------
    // Total mass [kg] = sum of layer masses.
    double calc_total_mass() const noexcept {
        double total = 0.0;
        for (const auto& layer : this->p_layers) { total += layer->get_mass(); }
        return total;
    }

    // Total internal radiogenic heating [W] at time_s. Only SolidLiquidLayers
    // produce radiogenic heating; other layer types contribute zero.
    double calc_internal_heating(double time_s) const noexcept {
        double total = 0.0;
        for (const auto& layer : this->p_layers) {
            const auto* sl = dynamic_cast<const c_SolidLiquidLayer*>(layer.get());
            if (sl != nullptr) {
                total += sl->calc_radiogenic_heating(time_s, sl->get_mass());
            }
        }
        return total;
    }

    // Returns true if every layer boundary is continuous (inner-to-outer) and the
    // innermost layer starts at radius 0.
    bool validate_layers() const noexcept {
        double prev_outer = 0.0;
        for (const auto& layer : this->p_layers) {
            const double inner = layer->get_radius_inner();
            const double tol   = d_LAYER_CONTINUITY_RTOL * (prev_outer > 1.0 ? prev_outer : 1.0);
            if (std::abs(inner - prev_outer) > tol) { return false; }
            prev_outer = layer->get_radius_outer();
        }
        return true;
    }

    // -----------------------------------------------------------------------
    // Binary I/O — world fields + layer count, then each layer's full record.
    // -----------------------------------------------------------------------
    void write_binary(std::ostream& out) const override {
        this->write_layered_binary(out, static_cast<uint32_t>(BinaryClassID::LayeredWorld));
    }

    void read_binary(std::istream& in, bool force = false) override {
        c_TidalPyBaseClass::read_binary(in, force);
        this->read_world_fields(in);
        uint64_t n_layers = 0;
        in.read(reinterpret_cast<char*>(&n_layers), sizeof(uint64_t));
        if (!in) {
            throw std::runtime_error("TidalPy: failed to read LayeredWorld binary data");
        }
        this->p_layers.clear();
        this->p_layers.reserve(n_layers);
        for (uint64_t i = 0; i < n_layers; ++i) {
            this->p_layers.push_back(c_layer_from_binary(in, force));
        }
    }

protected:
    // Shared writer so subclasses (e.g. c_GasGiantWorld) reuse the layout with
    // their own BinaryClassID.
    void write_layered_binary(std::ostream& out, uint32_t class_id) const {
        const uint64_t payload = this->world_payload_bytes() + sizeof(uint64_t);
        write_binary_header(out, class_id, payload);
        this->write_world_fields(out);
        const auto n_layers = static_cast<uint64_t>(this->p_layers.size());
        out.write(reinterpret_cast<const char*>(&n_layers), sizeof(uint64_t));
        if (!out) {
            throw std::runtime_error("TidalPy: failed to write layered-world binary data");
        }
        // Each layer writes its own complete record (recursively including its models).
        for (const auto& layer : this->p_layers) { layer->write_binary(out); }
    }

    std::vector<std::unique_ptr<c_BaseLayer>> p_layers;
};

} // namespace tidalpy
