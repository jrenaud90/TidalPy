# Obliquity Functions $F_{l,m,p}(I)$
_Note: The obliquity functions discussed in this section are often referred to as "inclination function" in the wider
literature. For the context of tides we want to use the angle between a planet/moon's spin axis and its orbital plane
around its primary host (the host would be Earth for our Moon or a star for an exoplanet). This angle is often called
usually described as the target planet/moon's obliquity angle. "Inclination" is an orbital parameter but is still
needed to find the relative obliquity of a host planet relative to its tidal target's orbit._

Obliquity functions are a key part of the Tidal Potential (See, e.g., the $F_{l,m,p}(I)$ in Eq. 1 of
[Kaula 1964](http://doi.wiley.com/10.1029/RG002i004p00661)). The potential is then used to find tidal strain, heating,
and orbit-spin evolution. It is critical to find the value of the obliquity function even if you assume your planet or
moon has 0 obliquity. Some $l,m,p$ tidal modes of the obliquity function are non-zero even if obliquity is.

TidalPy takes the approach of using pre-compiled functions that have already been analytically determined rather than
using the definition of $F_{l,m,p}(I)$. This greatly improves performance and provides a way to weed out tidal modes 
that have no impact on the potential (if $F_{l,m,p}(I) = 0$ for a particular $l,m,p$ then that mode does not
contribute to the potential).

To do this a 3rd-party script is run which builds a list of $F_{l,m,p}(I)$ for a given degree $l$ and truncation level.
The truncation level determines how many terms are included. Since many modes can still be active, it is best to pick
the lowest truncation level for your problem. The options currently available in TidalPy are:
- "general": No truncation takes place. This will be the most accurate solution but have the largest number of modes
    and therefore the slowest computation time.
- "4": The obliquity function is taylor expanded and truncated to only include $I^3$ and lower terms. A moderate amount
    of modes will be active. Accurate for moderate to small $I$.
- "2": The obliquity function is taylor expanded and truncated to only include $I$. A medium amount
    of modes will be active. Accurate for small $I$. This is the typical truncation found in the literature when
    non-zero obliquity is considered.
- "off": The obliquity function assumes that $I=0$. It is important to note that there will still be non-zero results
    for certain modes. The least amount of modes will be active. It is much better to use this truncation if you know
    $I=0$, it will greatly speed up calculations compared to using, for example, the general truncation and setting
    $I=0$.

Details about how to use the obliquity functions are below followed by a discussion of how they are analytically
determined.

## Using TidalPy's Obliquity Calculator
TidalPy provides obliquity function calculators for C++, Cython, and Python. Below we will go into how to use the
functions from Python. For users who would like to use them in C++ or Cython, we recommend looking into the files found
in `TidalPy.Tides_x.obliquity`. The "obliquity_common" files contain basic type information. "obliquity_driver" defines
the functions that will locate the correct C++ function to call based on user inputs for degree $l$ and truncation
level. All other files are the C++ (and Cython wrapped) obliquity functions for various truncations and degree $l$.

Using the Python Obliquity Calculator:
```python
from TidalPy.Tides_x.obliquity import obliquity_func

obliquity = np.radians(20.0)  # Obliquity must be passed in as radians.
results_l2 = obliquity_func(
    obliquity,
    degree_l=2        # First let's calculate for degree l=2
    truncation='gen'  # Let's say for l=2 we want the most accurate picture of the important modes.
)

results_l3 = obliquity_func(
    obliquity,
    degree_l=3      # Now for l=3
    truncation='2'  # Since the significance of modes inversely scales with degree l, perhaps we don't need as accurate a result for this l=3
)

# All obliquity results are stored in a pair of dict-like objects.
l2_results_by_lmp, l2_results_by_lm = results_l2

# `l2_results_by_lmp` is a dict-like object that stores obliquity results stored by a tuple of the mode integers l, m, p
print(f"Number of non-zero modes = {len(l2_results_by_lmp)}")
# The obliquity function _only returns non-zero modes_ all other modes might be present but result in $F_{lmp} = 0$ so are not returned by TidalPy's obliquity function.
print(f"Mode F_(l=2,m=0,p=0) = {l2_results_by_lmp[(2,0,0)]:0.3e}.)

# `l2_results_by_lm` is a dict-like object that stores further dict-like objects that contain the obliquity results.
# This can be useful if walking through unique (l, m) and want to find all results stored by p at that l,m
for (l, m), l2_results_by_p in l2_results_by_lm.items():
    for (p,), obliquity_result in l2_results_by_p:
        print(f"F_(l={l},m={m},p={p}) = {obliquity_result:0.3e}.")

# The format is the same for other degree l's and truncation levels. The only thing that will change is the number of
# non-zero modes and the result of the function.
```

## Determining Obliquity Functions Formulae
TidalPy using the following model for finding the obliquity functions. Note there is some debate on the nuances of some
steps in the process. We try to highlight those where applicable. If you see an error or have any questions please
leave a issue on TidalPy's GitHub repository.

The obliquity (or inclination) function is defined as:

$$F_{lmp}(I) = (-1)^{E(\frac{l-m+1}{2})} C_{lmp} S_{lmp}(I)$$

Where,

$$C_{lmp} = \frac{(l + m)!}{2^{l} p! (l-p)!}$$

and,

$$S_{lmp}(I) = \sum_{\lambda = \lambda_{min}}^{\lambda_{max}} (-1)^{\lambda} \binom{2l - 2p}{\lambda}\binom{2p}{l-m-\lambda}\left(\cos \frac{I}{2}\right)^{3l - m - 2p - 2\lambda}\left(\sin \frac{I}{2}\right)^{m -l + 2p + 2\lambda}$$.

The summation bounds are set by,

$$\lambda_{min} = \text{max}(0, l-m-2p)$$
$$\lambda_{max} = \text{min}(l-m, 2l - 2p)$$

The $E(x)$ function in the first equation is the "Entier" function or "integer part". It affectively acts as a floor 
function for positive $x$ and a ceiling function for negative $x$. This is often missing in a lot of publications,
it is usually not a big deal because for tidal heating we only need $F_{lmp}^{2}$ and squaring the function makes the
use of this function irrelevant. However, for the tidal potential we do not square $F$ so need to include it in our
definition here. See the conversation in
[Gooding & Wagner 2008](https://link.springer.com/article/10.1007/s10569-008-9145-6) and in the Appendix of
[Renaud et al 2021](https://iopscience.iop.org/article/10.3847/PSJ/abc0f3).

Unlike the eccentricity functions, the obliquity functions do not have infinite summations so can always be written
down exactly. This enables us to have the "general" truncation in TidalPy.
