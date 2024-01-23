# Rheology Functions and Classes
TidalPy provides an efficient suite of tools to calculate the complex shear modulus from a user-provided static shear, 
static viscosity, and static frequency. Note that "static" here does not mean "constant", instead it refers to the
purely real-valued versions of these parameters. All of them can change with time or radius.


## Defining a New Rheological Model
New rheologies can be implemented by creating a new sub class by cloning the TidalPy repository and adding a new
subclass to `\TidalPy\rheology\models.pyx`. You can copy one of the rheologies already present and use them as a
template to create your new model. After you are finished you will need to add the model to the
`find_rheology` function at the top of the same file as well as the `\TidalPy\rheology\models.pxd` header file 
so they can be imported by other TidalPy methods. You will then need to reinstall TidalPy so that it can compile
your new method (use `pip install -v .` from a terminal pointing to your TidalPy directory that you made the
modifications in). 

If you plan to push this new rheology to the main TidalPy Github: please consider adding test cases to 
`Tests\Test_Functions\test_rheology.py` so that TidalPy's CI system can check that it is performing correctly in 
future releases.