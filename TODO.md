
* Look into https://mybinder.org/ for when git goes unprivate.

* Add issue for logger not changing level if user changes tidalpy config after first load and calls reinit()
* Add looking into Cython or fortran for speed-up before 1.0 release to github issues
* Make known issues somewhere; mention numba.errors.TypingError and the flatten() workaround
* Once TidalPy is live on Github look into making a community on https://gitter.im/home/explore#
* Make a "build_from_planet" function to take an already built planet and then provide some global changes that builds a new planet with, for example, a larger mass/radius
    * Make it so if only mass, or only radius, is provided an approx (M-R) of the other is used?
* add readme on how to get jupyter widgets working: https://towardsdatascience.com/interactive-controls-for-jupyter-notebooks-f5c94829aee6
## On Future Python or 3rd Party Package Updates
* If python's builtin `json` package updates to support json5 then it should be used over the currently supported 3rd party `json5` package.
* Numba currently does not support array slicing with 2D index arrays.
    * Currently does not work if `x` is 2D:
    ```
        res = x**t
        res[res<2.] = 2.
    ```
    * Workaround that is njit-safe
    ```
        res = x**t
        s = res.shape
        res_f = res.flatten()
        res_f[res_f<2.] = 2.
        res_f = np.reshape(res_f, s)
    ``` 
  
* Before major release:
    * Change default logging level to info