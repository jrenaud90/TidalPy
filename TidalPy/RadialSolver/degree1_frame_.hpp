#pragma once
/* Reference frames of the degree-1 load Love numbers (Blewitt 2003).
 *
 * At degree 1 a rigid translation of the body meets every surface condition of a load, so the load Love numbers are
 * defined only once the origin is chosen. The radial solver computes them in CE, the center of mass of the solid
 * body (k' = 0), and shifts them to another frame afterwards. A shift of the origin along the load moment adds the
 * same amount to h', l', and k' (Blewitt 2003 Eq. 17): X_frame = X_CE - alpha for X in h', l', and 1 + k', with
 *     CE  center of mass of the solid body                      alpha = 0
 *     CM  center of mass of the body and the load (1 + k' = 0)  alpha = 1
 *     CF  center of surface figure (h' + 2 l' = 0)              alpha = (h' + 2 l') / 3
 *     CL  center of lateral figure (l' = 0)                     alpha = l'
 *     CH  center of height figure (h' = 0)                      alpha = h'
 * with h' and l' in CE (Blewitt 2003 Eqs. 18 to 25; Martens 2016 Eqs. 4.202 and 4.203). In the radial functions the
 * shift is the translation y1 = y3 = c, y5 = c g(r), with c = -alpha / g at the surface.
 */

#include <complex>
#include <limits>


enum class c_Degree1Frame : int
{
    CE = 0,
    CM = 1,
    CF = 2,
    CL = 3,
    CH = 4
};

inline constexpr int C_NUM_DEGREE1_FRAMES = 5;


// Whether the frame's offset reads l', which a static liquid surface does not define.
inline constexpr bool c_degree1_frame_needs_l(int frame) noexcept
{
    return (frame == static_cast<int>(c_Degree1Frame::CF)) || (frame == static_cast<int>(c_Degree1Frame::CL));
}


// Blewitt's alpha: the amount h', l', and 1 + k' in CE drop by in the frame. NaN for an unknown frame.
inline std::complex<double> c_degree1_frame_offset(
        int frame,
        const std::complex<double>& h_ce,
        const std::complex<double>& l_ce) noexcept
{
    switch (static_cast<c_Degree1Frame>(frame))
    {
        case c_Degree1Frame::CE: return std::complex<double>(0.0, 0.0);
        case c_Degree1Frame::CM: return std::complex<double>(1.0, 0.0);
        case c_Degree1Frame::CF: return (h_ce + 2.0 * l_ce) / 3.0;
        case c_Degree1Frame::CL: return l_ce;
        case c_Degree1Frame::CH: return h_ce;
    }
    return std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
}
