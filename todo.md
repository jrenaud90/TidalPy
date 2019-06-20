## Planets
* Have the dill loader check the TidalPy version that a planet was saved under and either raise an error or a warning if user tries to load an incompatible version
* Make it so that a layer can have just 1 rheology and that different rheologies can collapse onto one another. e.g., Mantle = Maxwell, Crust = Andrade -> Global love = sum of the two different rheologies (instead of storing them separately in a dictionary)
* Allow planets to load with a default eccentricity and inclination

## Layers
* As of 0.1.0 Stefan number is not calculated using melt fraction (just a configuration constant)

## Thermal
* Move the viscosity calculation (pre-melt) into rheology module (currently in thermal)?

## Rheology
* Put a check to make sure that alpha stays between 0 and 1 (the factorial code will not work outside that range)

## Code
* Change the model search to look for "Parameters\n----------" in doc string and then each parameter will have a flag: "**input" "**live" or similar
* Make the setters in orbit class accept a list for the planet ref so that multiple planets can be changed at once. 