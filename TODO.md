* Make collaboration guide with notes on TODO, FIXME, OPT
* Make known issues somewhere; mention numba.errors.TypingError and the flatten() workaround
* Once TidalPy is live on Github look into making a community on https://gitter.im/home/explore#
* Make a "build_from_planet" function to take an already built planet and then provide some global changes that builds a new planet with, for example, a larger mass/radius
    * Make it so if only mass, or only radius, is provided an approx (M-R) of the other is used?
* Make a custom vectorize function that uses numba.vectorize if numba enables otherwise uses np.vectorize.
  This will be a decorator with arguments so that type information can be passed to the numba vectorize - it looks like this signiture is not required but could be provide a speed up.
* Go through and turn back on many of the commented `#njit` commands.
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