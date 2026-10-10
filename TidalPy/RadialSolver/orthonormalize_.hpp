#pragma once
/* Re-orthonormalization of a layer's independent shooting solutions (Godunov 1961; Conte 1966).
 *
 * The solutions a layer starts from grow at different rates, so they turn toward the fastest growing one and the
 * combinations the boundary conditions need cancel. The surface solve then loses digits in proportion to how
 * dependent they have become. The shooting method integrates a layer's solutions as one system with a terminal CyRK
 * event on their normalized Gram determinant (c_solution_independence_event). Where the integration has reduced it by
 * the factor [numerical] minimum_solution_independence, the integration stops, the solutions are replaced by an
 * orthonormal basis of the same span, Y = Q R (c_orthonormalize_solutions), and a new integration continues from it.
 * The ODE is linear, so the solutions before the restart are those after it times R, Y_before(r) = Y_after(r) R, at
 * every radius, and the collapse constants of the earlier segment are c_before = R^-1 c_after
 * (c_constants_before_basis_change). A dynamic liquid's first segment holds another basis change, the split of its
 * starting solutions into one with P exactly 0 and one that carries P (c_split_pressure_solutions), which need not
 * be triangular.
 *
 * Both the independence and the new basis are taken in the ys divided by their characteristic sizes, TidalPy's
 * non-dimensional units (times the weights of c_layer_y_weights), so a solve in SI behaves as a non-dimensional one. In
 * SI a stress is about 1e10 times a displacement, and the raw ys would look dependent as soon as they were
 * orthonormal.
 *
 * Solutions are stored as CyRK holds them: num_solutions blocks of num_ys complex ys, two reals each.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>

#include <Eigen/Dense>

#include "rs_constants_.hpp"        // C_MAX_NUM_SOL, C_MAX_NUM_Y, c_BasisChange
#include "derivatives/odes_.hpp"    // c_RadialSolverArgs


/// Normalized Gram determinant of the scaled solutions w_ik = y_ik y_weights[k], det(G) / prod(G_ii) with
/// G_ij = sum_k conj(w_ik) w_jk: 1 when they are orthogonal, 0 when they are dependent (or one vanishes), and 1 for a
/// single solution. Each solution is divided by its largest scaled component first, so the products cannot overflow.
inline double c_solution_independence(
        const double* y_ptr,
        size_t num_solutions,
        size_t num_ys,
        const double* y_weights_ptr) noexcept
{
    if (num_solutions < 2)
    {
        return 1.0;
    }
    std::complex<double> scaled[C_MAX_NUM_SOL][C_MAX_NUM_Y];
    for (size_t solution_i = 0; solution_i < num_solutions; ++solution_i)
    {
        const double* solution_ptr = y_ptr + 2 * num_ys * solution_i;
        double largest = 0.0;
        for (size_t y_i = 0; y_i < num_ys; ++y_i)
        {
            scaled[solution_i][y_i] =
                std::complex<double>(solution_ptr[2 * y_i], solution_ptr[2 * y_i + 1]) * y_weights_ptr[y_i];
            largest = std::max(largest, std::abs(scaled[solution_i][y_i]));
        }
        if (!(largest > 0.0) || !std::isfinite(largest))
        {
            return 0.0;
        }
        const double inverse_largest = 1.0 / largest;
        for (size_t y_i = 0; y_i < num_ys; ++y_i)
        {
            scaled[solution_i][y_i] *= inverse_largest;
        }
    }
    std::complex<double> gram[C_MAX_NUM_SOL][C_MAX_NUM_SOL];
    for (size_t row_i = 0; row_i < num_solutions; ++row_i)
    {
        for (size_t column_i = row_i; column_i < num_solutions; ++column_i)
        {
            std::complex<double> product(0.0, 0.0);
            for (size_t y_i = 0; y_i < num_ys; ++y_i)
            {
                product += std::conj(scaled[row_i][y_i]) * scaled[column_i][y_i];
            }
            gram[row_i][column_i] = product;
        }
    }
    const double g00 = gram[0][0].real();
    const double g11 = gram[1][1].real();
    if (num_solutions == 2)
    {
        return 1.0 - std::norm(gram[0][1]) / (g00 * g11);
    }
    const double g22 = gram[2][2].real();
    // The determinant of a 3 x 3 Hermitian matrix, which is real.
    const double determinant = g00 * g11 * g22
        + 2.0 * (gram[0][1] * gram[1][2] * std::conj(gram[0][2])).real()
        - g00 * std::norm(gram[1][2]) - g11 * std::norm(gram[0][2]) - g22 * std::norm(gram[0][1]);
    return determinant / (g00 * g11 * g22);
}


/// CyRK EventFunc: the solutions' independence less the args' floor, which falls through zero where the shooting
/// method re-orthonormalizes them.
inline double c_solution_independence_event(double radius, double* y_ptr, char* args_ptr) noexcept
{
    (void)radius;
    const c_RadialSolverArgs* rs_args_ptr = reinterpret_cast<const c_RadialSolverArgs*>(args_ptr);
    const double independence = c_solution_independence(
        y_ptr, rs_args_ptr->num_solutions, rs_args_ptr->num_ys, rs_args_ptr->y_weights.data());
    return independence - rs_args_ptr->independence_floor;
}


/// Replace the solutions with an orthonormal basis of their span in the scaled ys, in place. With W the diagonal of
/// y_weights, the Householder QR of W Y = Q R gives the new solutions W^-1 Q. basis_change receives R (upper
/// triangular), so that the solutions given are the new ones times R.
inline void c_orthonormalize_solutions(
        double* y_ptr,
        size_t num_solutions,
        size_t num_ys,
        const double* y_weights_ptr,
        c_BasisChange& basis_change) noexcept
{
    using c_SolutionMatrix = Eigen::Matrix<
        std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, C_MAX_NUM_Y, C_MAX_NUM_SOL>;
    c_SolutionMatrix scaled(num_ys, num_solutions);
    for (size_t solution_i = 0; solution_i < num_solutions; ++solution_i)
    {
        for (size_t y_i = 0; y_i < num_ys; ++y_i)
        {
            const size_t value_i = 2 * (num_ys * solution_i + y_i);
            scaled(y_i, solution_i) = std::complex<double>(y_ptr[value_i], y_ptr[value_i + 1]) * y_weights_ptr[y_i];
        }
    }
    const Eigen::HouseholderQR<c_SolutionMatrix> qr(scaled);
    const c_SolutionMatrix orthonormal = qr.householderQ() * c_SolutionMatrix::Identity(num_ys, num_solutions);

    basis_change.fill(std::complex<double>(0.0, 0.0));
    for (size_t row_i = 0; row_i < num_solutions; ++row_i)
    {
        for (size_t column_i = row_i; column_i < num_solutions; ++column_i)
        {
            basis_change[row_i * C_MAX_NUM_SOL + column_i] = qr.matrixQR()(row_i, column_i);
        }
    }
    for (size_t solution_i = 0; solution_i < num_solutions; ++solution_i)
    {
        for (size_t y_i = 0; y_i < num_ys; ++y_i)
        {
            const size_t value_i = 2 * (num_ys * solution_i + y_i);
            const std::complex<double> value = orthonormal(y_i, solution_i) / y_weights_ptr[y_i];
            y_ptr[value_i]     = value.real();
            y_ptr[value_i + 1] = value.imag();
        }
    }
}


/// A dynamic liquid's two starting solutions recombined, in place, into one whose P is exactly 0 and one that carries
/// the layer's P: with s_j the solution of larger P, s_j scaled to a largest weighted component of 1, and
/// s_k - (P_k / P_j) s_j with its P set to 0. A physical solution's P is of order omega^2 (c_pressure_unit), and the
/// integration carries it to y1 through y3 = -P / (rho omega^2 r), so the combination of the starting solutions that
/// cancels P is formed here, before its roundoff is divided by omega^2. basis_change receives T, so that the solutions
/// given are the new ones times T.
inline void c_split_pressure_solutions(
        double* y_ptr,
        size_t num_ys,
        const double* y_weights_ptr,
        c_BasisChange& basis_change) noexcept
{
    basis_change = c_identity_basis_change();
    const auto value = [&](size_t solution_i, size_t y_i)
    {
        const size_t value_i = 2 * (num_ys * solution_i + y_i);
        return std::complex<double>(y_ptr[value_i], y_ptr[value_i + 1]);
    };
    const auto store = [&](size_t solution_i, size_t y_i, const std::complex<double>& new_value)
    {
        const size_t value_i = 2 * (num_ys * solution_i + y_i);
        y_ptr[value_i]     = new_value.real();
        y_ptr[value_i + 1] = new_value.imag();
    };
    const size_t pressure_i = C_DYNAMIC_LIQUID_LAYOUT.slot_of_full_y(1);   // P is stored in y2's slot
    const size_t carrier_i  = (std::abs(value(1, pressure_i)) > std::abs(value(0, pressure_i))) ? 1 : 0;
    const size_t other_i    = 1 - carrier_i;
    const std::complex<double> carrier_pressure = value(carrier_i, pressure_i);
    if (!(std::abs(carrier_pressure) > 0.0))
    {
        return;
    }
    const std::complex<double> ratio = value(other_i, pressure_i) / carrier_pressure;
    double largest = 0.0;
    for (size_t y_i = 0; y_i < num_ys; ++y_i)
    {
        store(other_i, y_i, (y_i == pressure_i) ? std::complex<double>(0.0, 0.0) :
            value(other_i, y_i) - ratio * value(carrier_i, y_i));
        largest = std::max(largest, std::abs(value(carrier_i, y_i)) * y_weights_ptr[y_i]);
    }
    for (size_t y_i = 0; y_i < num_ys; ++y_i)
    {
        store(carrier_i, y_i, value(carrier_i, y_i) / largest);
    }
    basis_change[carrier_i * C_MAX_NUM_SOL + carrier_i] = largest;
    basis_change[carrier_i * C_MAX_NUM_SOL + other_i]   = ratio * largest;
}


/// The basis change of two in turn: solutions Y = Y' first and Y' = Y'' second give Y = Y'' (second first).
inline c_BasisChange c_compose_basis_changes(
        const c_BasisChange& first,
        const c_BasisChange& second,
        size_t num_solutions) noexcept
{
    c_BasisChange product{};
    for (size_t row_i = 0; row_i < num_solutions; ++row_i)
    {
        for (size_t column_i = 0; column_i < num_solutions; ++column_i)
        {
            std::complex<double> sum(0.0, 0.0);
            for (size_t inner_i = 0; inner_i < num_solutions; ++inner_i)
            {
                sum += second[row_i * C_MAX_NUM_SOL + inner_i] * first[inner_i * C_MAX_NUM_SOL + column_i];
            }
            product[row_i * C_MAX_NUM_SOL + column_i] = sum;
        }
    }
    return product;
}


/// The constants of the solutions before a basis change from those after it, c_before = B^-1 c_after (Gaussian
/// elimination with partial pivoting; B is the upper-triangular R of a re-orthonormalization, or a dynamic liquid's
/// starting split). The two arrays may be the same.
inline void c_constants_before_basis_change(
        const c_BasisChange& basis_change,
        const std::complex<double>* constants_after_ptr,
        std::complex<double>* constants_before_ptr,
        size_t num_solutions) noexcept
{
    std::complex<double> matrix[C_MAX_NUM_SOL][C_MAX_NUM_SOL];
    std::complex<double> rhs[C_MAX_NUM_SOL];
    for (size_t row_i = 0; row_i < num_solutions; ++row_i)
    {
        for (size_t column_i = 0; column_i < num_solutions; ++column_i)
        {
            matrix[row_i][column_i] = basis_change[row_i * C_MAX_NUM_SOL + column_i];
        }
        rhs[row_i] = constants_after_ptr[row_i];
    }
    for (size_t pivot_i = 0; pivot_i < num_solutions; ++pivot_i)
    {
        size_t best_i = pivot_i;
        for (size_t row_i = pivot_i + 1; row_i < num_solutions; ++row_i)
        {
            if (std::abs(matrix[row_i][pivot_i]) > std::abs(matrix[best_i][pivot_i])) { best_i = row_i; }
        }
        if (best_i != pivot_i)
        {
            for (size_t column_i = 0; column_i < num_solutions; ++column_i)
            {
                std::swap(matrix[pivot_i][column_i], matrix[best_i][column_i]);
            }
            std::swap(rhs[pivot_i], rhs[best_i]);
        }
        for (size_t row_i = pivot_i + 1; row_i < num_solutions; ++row_i)
        {
            const std::complex<double> factor = matrix[row_i][pivot_i] / matrix[pivot_i][pivot_i];
            if (factor == std::complex<double>(0.0, 0.0)) { continue; }
            for (size_t column_i = pivot_i; column_i < num_solutions; ++column_i)
            {
                matrix[row_i][column_i] -= factor * matrix[pivot_i][column_i];
            }
            rhs[row_i] -= factor * rhs[pivot_i];
        }
    }
    for (size_t row_i = num_solutions; row_i-- > 0;)
    {
        std::complex<double> sum = rhs[row_i];
        for (size_t column_i = row_i + 1; column_i < num_solutions; ++column_i)
        {
            sum -= matrix[row_i][column_i] * constants_before_ptr[column_i];
        }
        constants_before_ptr[row_i] = sum / matrix[row_i][row_i];
    }
}
