"""Generate the angular Gram table for the analytic colatitude collapse of the secular 3D heating.

The secular 3D tidal volumetric heating of a single (l, m) mode factorizes as
    h_mode(r, theta) / amplitude^2 = sum_k w_k Im( sigma_tilde_k(r, theta) conj(eps_tilde_k(r, theta)) )
where each reduced strain/stress component is (complex radial coefficient) x (real angular function of
theta). The six angular functions form a bounded basis (the 1/sin(theta) poles cancel in the phi-phi and
theta-phi combinations):

    f1 = P_lm(cos theta)
    f2 = dP/dtheta
    f3 = d2P/dtheta2
    f4 = P / sin theta
    f5 = -m^2 P / sin^2 theta + cot(theta) dP/dtheta          (the phi-phi angular operator)
    f6 = ( dP/dtheta - cot(theta) P ) / sin theta             (the theta-phi angular operator)

Collapsing the colatitude integral analytically then needs only the symmetric Gram matrix
    G_ij(l, m) = integral_0^pi f_i(theta) f_j(theta) sin(theta) dtheta
so that
    integral_0^pi h_mode sin(theta) dtheta / amplitude^2
        = sum_k w_k sum_ij Im( Ssigma_ki conj(Seps_kj) ) G_ij(l, m).

This script computes G_ij with sympy (exact P_lm, Condon-Shortley phase, matching Utilities_x/legendre)
and mpmath high-precision quadrature, and writes them as double constants to angular_gram_.hpp. Run:

    python -m TidalPy.Tides_x.multilayer.codegen.gen_angular_gram   (or run this file directly)

then rebuild. The generated table is validated independently in
Tests/Test_Tides_x/Test_Multilayer/test_angular_gram_01.py.
"""
import os

import mpmath as mp
import sympy as sp

MIN_L = 2
MAX_L = 10
BASIS_SIZE = 6


def basis_functions(degree_l, order_m):
    """Return the six bounded angular basis functions f1..f6 as sympy expressions in x = cos(theta).

    Working in x rather than theta keeps each integrand a smooth (polynomial x sqrt(1-x^2)) function, so
    the mpmath tanh-sinh quadrature is accurate. (In theta the high-degree derivatives oscillate rapidly
    and adaptive quadrature badly under-resolves them.) The theta-derivatives are expressed via
    dP/dtheta = -sin dP/dx and d2P/dtheta2 = sin^2 d2P/dx2 - cos dP/dx, with sin = sqrt(1-x^2), cos = x.
    """
    x = sp.symbols("x", real=True)
    sin_t = sp.sqrt(1 - x ** 2)
    cot_t = x / sin_t
    legendre = sp.assoc_legendre(degree_l, order_m, x)
    d_legendre_dx = sp.diff(legendre, x)
    d_legendre = -sin_t * d_legendre_dx                                   # dP/dtheta
    d2_legendre = (1 - x ** 2) * sp.diff(legendre, x, 2) - x * d_legendre_dx  # d2P/dtheta2
    funcs = [
        legendre,
        d_legendre,
        d2_legendre,
        legendre / sin_t,
        -order_m ** 2 * legendre / sin_t ** 2 + cot_t * d_legendre,
        (d_legendre - cot_t * legendre) / sin_t,
    ]
    # Simplify so the removable 1/sin poles are cancelled before lambdifying (f5, f6 are bounded).
    return x, [sp.simplify(func) for func in funcs]


# For m = 0 the f4 = P/sin and f6 = (dP - cot*P)/sin basis functions are UNBOUNDED (P_l0 does not
# vanish at the poles), so their Gram integrals diverge. They enter the strain only through the factor
# i*m (eps_rphi ~ i m f4, eps_thphi ~ i m f6), which is exactly zero for m = 0, so those entries are
# physically unused; we set them to 0 rather than store the divergent (garbage) quadrature. For m >= 1
# every basis function is bounded and every entry is a genuine finite integral.
_ZERO_FOR_M0 = {3, 5}  # 0-based indices of f4 and f6


def _quad(func_i, func_j, dps):
    """int_{-1}^{1} f_i(x) f_j(x) dx (= int_0^pi f_i f_j sin(theta) dtheta) by tanh-sinh quadrature.

    Split at 0 so the (integrable sqrt(1-x^2)) endpoints at +-1 are each an interval endpoint.
    """
    mp.mp.dps = dps
    return mp.quad(lambda xx: func_i(xx) * func_j(xx), [-1, 0, 1])


def gram_matrix(degree_l, order_m):
    """The symmetric 6x6 Gram matrix G_ij = int_0^pi f_i f_j sin(theta) dtheta.

    Uses adaptive tanh-sinh quadrature in x = cos(theta) at two precisions, with a convergence self-check
    so a stored value can never be a silently-unconverged quadrature.
    """
    x, funcs = basis_functions(degree_l, order_m)
    lambdified = [sp.lambdify(x, func, "mpmath") for func in funcs]
    gram = [[0.0] * BASIS_SIZE for _ in range(BASIS_SIZE)]
    for i in range(BASIS_SIZE):
        for j in range(i, BASIS_SIZE):
            if order_m == 0 and (i in _ZERO_FOR_M0 or j in _ZERO_FOR_M0):
                continue  # f4/f6 diverge for m = 0 but are multiplied by i*m = 0 downstream (unused)
            coarse = _quad(lambdified[i], lambdified[j], 30)
            fine = _quad(lambdified[i], lambdified[j], 50)
            scale = max(abs(fine), abs(coarse), mp.mpf(1))
            if abs(fine - coarse) / scale > mp.mpf("1e-12"):
                raise RuntimeError(
                    f"Gram integral not converged at l={degree_l}, m={order_m}, "
                    f"(i,j)=({i},{j}): {coarse!r} vs {fine!r}")
            gram[i][j] = gram[j][i] = float(fine)
    return gram


def upper_triangle(gram):
    """The 21 upper-triangle (i <= j) entries row-major."""
    return [gram[i][j] for i in range(BASIS_SIZE) for j in range(i, BASIS_SIZE)]


def generate(header_path):
    rows = []
    labels = []
    for degree_l in range(MIN_L, MAX_L + 1):
        for order_m in range(0, degree_l + 1):
            rows.append(upper_triangle(gram_matrix(degree_l, order_m)))
            labels.append((degree_l, order_m))
            print(f"  computed Gram for l={degree_l}, m={order_m}")

    num_pairs = len(rows)
    lines = []
    lines.append("#pragma once")
    lines.append("// GENERATED by codegen/gen_angular_gram.py - do not edit by hand.")
    lines.append("//")
    lines.append("// Angular Gram table for the analytic colatitude collapse of the secular 3D tidal heating.")
    lines.append("// G_ij(l, m) = int_0^pi f_i(theta) f_j(theta) sin(theta) dtheta for the bounded 6-function basis")
    lines.append("//   f1 = P_lm, f2 = dP/dtheta, f3 = d2P/dtheta2, f4 = P/sin,")
    lines.append("//   f5 = -m^2 P/sin^2 + cot*dP,  f6 = (dP - cot*P)/sin.")
    lines.append("// Symmetric 6x6; the 21 upper-triangle entries (i <= j, row-major) are stored per (l, m),")
    lines.append("// for l = 2..10 and m = 0..l. See gen_angular_gram.py for the derivation.")
    lines.append("#include <cstddef>")
    lines.append("")
    lines.append("namespace tidalpy {")
    lines.append("namespace tides {")
    lines.append("")
    lines.append(f"inline constexpr int C_ANGULAR_GRAM_MIN_L = {MIN_L};")
    lines.append(f"inline constexpr int C_ANGULAR_GRAM_MAX_L = {MAX_L};")
    lines.append(f"inline constexpr int C_ANGULAR_GRAM_NUM_PAIRS = {num_pairs};")
    lines.append("inline constexpr int C_ANGULAR_GRAM_BASIS = 6;")
    lines.append("inline constexpr int C_ANGULAR_GRAM_UPPER = 21;")
    lines.append("")
    lines.append("// Row index for (l, m): sum_{l'=2}^{l-1}(l'+1) + m.")
    lines.append("inline int c_angular_gram_row(int degree_l, int order_m) {")
    lines.append("    int row = 0;")
    lines.append("    for (int ll = C_ANGULAR_GRAM_MIN_L; ll < degree_l; ++ll) { row += (ll + 1); }")
    lines.append("    return row + order_m;")
    lines.append("}")
    lines.append("")
    lines.append("// Upper-triangle (i <= j, row-major) Gram entries per (l, m).")
    lines.append("inline constexpr double C_ANGULAR_GRAM_TABLE[C_ANGULAR_GRAM_NUM_PAIRS][C_ANGULAR_GRAM_UPPER] = {")
    for (degree_l, order_m), row in zip(labels, rows):
        entries = ", ".join(f"{value: .17e}" for value in row)
        lines.append(f"    {{ {entries} }},  // l={degree_l}, m={order_m}")
    lines.append("};")
    lines.append("")
    lines.append("// Fill a full symmetric 6x6 Gram matrix for (l, m). Returns false if (l, m) is out of range.")
    lines.append("inline bool c_angular_gram(int degree_l, int order_m, double gram[6][6]) {")
    lines.append("    if (degree_l < C_ANGULAR_GRAM_MIN_L || degree_l > C_ANGULAR_GRAM_MAX_L")
    lines.append("        || order_m < 0 || order_m > degree_l) {")
    lines.append("        return false;")
    lines.append("    }")
    lines.append("    const double* upper = C_ANGULAR_GRAM_TABLE[c_angular_gram_row(degree_l, order_m)];")
    lines.append("    int idx = 0;")
    lines.append("    for (int i = 0; i < 6; ++i) {")
    lines.append("        for (int j = i; j < 6; ++j) {")
    lines.append("            gram[i][j] = gram[j][i] = upper[idx++];")
    lines.append("        }")
    lines.append("    }")
    lines.append("    return true;")
    lines.append("}")
    lines.append("")
    lines.append("}  // namespace tides")
    lines.append("}  // namespace tidalpy")
    lines.append("")

    with open(header_path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))
    print(f"wrote {header_path} ({num_pairs} (l, m) pairs)")


if __name__ == "__main__":
    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "angular_gram_.hpp")
    generate(out)
