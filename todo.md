## Planets
* Have the dill loader check the TidalPy version that a planet was saved under and either raise an error or a warning if user tries to load an incompatible version
* Make it so that a layer can have just 1 rheology and that different rheologies can collapse onto one another. e.g., Mantle = Maxwell, Crust = Andrade -> Global love = sum of the two different rheologies (instead of storing them separately in a dictionary)
* Allow planets to load with a default eccentricity and inclination

## Layers
* As of 0.1.0 Stefan number is not calculated using melt fraction (just a configuration constant)

## Thermal
* Move the viscosity calculation (pre-melt) into rheology module (currently in thermal)? 