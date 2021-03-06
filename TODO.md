* rheology.viscosity.viscosity_models.py still has `_array` functions.
* rheology.partial_melt.melting_models.py still has `_array` functions.
* rheology.partial_melt.partialmelt.py still has `_array` functions.
* A multilayer model has now been implemented, but it is not currently part of the OOP scheme.
* Probably worth refactoring `burnman_interface` into the `utilities` module to match other 3rd party packages.
* Add issue for logger not changing level if user changes tidalpy config after first load and calls reinit()
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