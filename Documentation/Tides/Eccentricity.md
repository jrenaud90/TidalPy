# Obliquity Functions $G_{l,p,q}(e)$
Similar to the Obliquity Functions, Eccentricity Functions are required to calculate the Tidal Potential (See, e.g.,
the $G_{l,m,p}(e)$ in Eq. 1 of [Kaula 1964](http://doi.wiley.com/10.1029/RG002i004p00661)). Unlike the obliquity
functions, the eccentricity functions are defined by an infinite sum ($q \in \{-\inf, \inf\}$) and can not
be written down exactly. This requires us to pick a truncation level which, as long as $e<1$, allows us to balance
errors with calculation speeds.

Like with the Obliquity Functions, TidalPy takes the approach of pre-calculating the terms in the Eccentricity
Functions as solving for them directly is very computationally expensive. A 3rd party script generates the terms for
a given degree l and truncation level. While we can go to arbitrary truncation level, we find that certain steps are
most important and are the ones that are currently provided with TidalPy. The number $n$ in brackets indicates
the number of tidal modes when $l=2$ that truncation activates giving you a sense for the increase in computation. The
higher the degree $l$ the more modes will be activated.
- Truncation=1 ($n=3$), or truncating to terms with powers $e^{1}$ closely matches the traditional tidal heating formula. The
    tidal potential is squared when finding heating, so truncating to powers <2 will result in $e^{2}$ terms 
    in the heating equation. Keep in mind that this is a bit more accurate since we keep exact terms like 
    $(1- e^2)^(-1.5)$ which does not get truncated. While this is more accurate, it does cause solutions to not
    exactly match the traditional $e^2$ formulation.
- Truncation=2 ($n=9$), truncating potential terms to $e^{2}$ (heating to $e^{4}$). Popular truncation.
- Truncation=3 ($n=13$), truncating potential terms to $e^{3}$ (heating to $e^{6}$). Popular truncation.
- Truncation=4 ($n=19$), truncating potential terms to $e^{4}$ (heating to $e^{8}$).
- Truncation=5 ($n=25$), truncating potential terms to $e^{5}$ (heating to $e^{10}$).
- Truncation=10 ($n=55$), truncating potential terms to $e^{10}$ (heating to $e^{20}$). Used in [Renaud et al 2021](https://iopscience.iop.org/article/10.3847/PSJ/abc0f3).
- Truncation=15 ($n=85$), truncating potential terms to $e^{15}$ (heating to $e^{30}$). I am curious so I put it in.
- Truncation=20 ($n=115$), truncating potential terms to $e^{20}$ (heating to $e^{40}$). If you like the smell of burnt silicon.

Some Notes:
- Some of these really large truncations are likely to have numerical errors. You should take the results they produce
    with a grain of salt. Perform some sensitivity testing by increasing the value of $e$ into the region where these
    truncations begin to matter. See if the results look stable.
- While some of these truncations seem crazy, Renaud et al. 2021 found that you need to get into the $e^{10}$ range
    once your eccentricity gets larger than 0.5. For example, tidal heating at truncation=5 is up to 2 orders of
    magnitude lower at $e=0.8$ compared to truncation=10.
- Keep in mind that increasing the truncation level activates new tidal modes. New tidal modes may unlock new unique
    forcing frequencies which will lead to new Love number calculations. So doubling the number of modes does not
    simply double the computation time, it is exponential. This is even worse at degree $l>2$ as higher degrees
    also unlock new modes.

Details about how to use the Eccentricity Functions are below followed by a discussion of how they are analytically
determined.

## Using TidalPy's Eccentricity Functions
TidalPy provides Eccentricity Function calculators for C++, Cython, and Python. Below we will go into how to use the
functions from Python. For users who would like to use them in C++ or Cython, we recommend looking into the files found
in `TidalPy.Tides_x.eccentricity`. The "eccentricity_common" files contain basic type information.
"eccentricity_driver" defines the functions that will locate the correct C++ function to call based on user inputs
for degree $l$ and truncation level. All other files are the C++ (and Cython wrapped) eccentricity functions for various
truncations and degree $l$.

Using the Python Eccentricity Calculator:
```python
from TidalPy.Tides_x.eccentricity import eccentricity_func

eccentricity = 0.1  # eccentricity must be greater than or equal to 0 and less than 1.
results_l2 = eccentricity_func(
    eccentricity,
    degree_l=2    # First let's calculate for degree l=2
    truncation=4  # Let's say for l=2 we want a pretty accurate picture of the important modes.
)

# As of TidalPy v0.7.x - degrees l = 2, ..., 10 are supported.
results_l3 = eccentricity_func(
    eccentricity,
    degree_l=3      # Now for l=3
    truncation=2  # Since the significance of modes inversely scales with degree l, perhaps we don't need as accurate a result for this l=3
)

# All eccentricity results are stored in a pair of dict-like objects.
l2_results_by_lpq, l2_results_by_lp = results_l2

# `l2_results_by_lpq` is a dict-like object that stores eccentricity results stored by a tuple of the mode integers l, p, q
print(f"Number of non-zero modes = {len(l2_results_by_lpq)}")
# The eccentricity function _only returns non-zero modes_ all other modes might be present but result in $G_{lmp} = 0$ so are not returned by TidalPy's eccentricity function.
print(f"Mode F_(l=2,m=0,p=0) = {l2_results_by_lpq[(2,0,0)]:0.3e}.)

# `l2_results_by_lp` is a dict-like object that stores further dict-like objects that contain the eccentricity results.
# This can be useful if walking through unique (l, m) and want to find all results stored by q at that l,p
for (l, p), l2_results_by_q in l2_results_by_lp.items():
    for (q,), eccentricity_result in l2_results_by_q:
        print(f"F_(l={l},p={p},q={q}) = {eccentricity_result:0.3e}.")

# The format is the same for other degree l's and truncation levels. The only thing that will change is the number of
# non-zero modes and the result of the function.
```

## Analytical Formulation of the Eccentricity Functions
The eccentricity functions have a long history in astrophysics and many different ways to calculate them. TidalPy
follows the approach discussed in Appendix C.1 of
[Renaud et al. 2021](https://iopscience.iop.org/article/10.3847/PSJ/abc0f3) the interested reader should also review
the work of [Hughes 1981](https://ui.adsabs.harvard.edu/abs/1981CeMec..25..101H/abstract),
[Laskar and Boue 2010](https://www.aanda.org/articles/aa/full_html/2010/14/aa14496-10/aa14496-10.html), and
[Veras et al. 2019](https://ui.adsabs.harvard.edu/abs/2019MNRAS.486.3831V/abstract).

In short, the eccentricity functions can be related to the "Hansen coefficients" via,
$$G_{lpq}(e) = H_{l-2p+q}^{-l-1,l-2p}(e)$$.

The Hansen coefficients are calculated in a number of different regimes depending on the values of $l,p,q$ (details
in Appendix C of Renaud et al. 2021). There are many techniques in the literature to improve the computational
efficiency of finding these coefficients. However, we choose mathematical accuracy over computational speed since the
functional form of the eccentricity functions are found ahead of time and pre-compiled.
