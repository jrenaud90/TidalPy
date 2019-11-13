* Make a "build_from_planet" function to take an already built planet and then provide some global changes that builds a new planet with, for example, a larger mass/radius
    * Make it so if only mass, or only radius, is provided an approx (M-R) of the other is used?
* Make a custom vectorize function that uses numba.vectorize if numba enables otherwise uses np.vectorize.
  This will be a decorator with arguments so that type information can be passed to the numba vectorize - it looks like this signiture is not required but could be provide a speed up.