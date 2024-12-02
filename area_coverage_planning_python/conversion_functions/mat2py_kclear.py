# This code implements the functions that call the SPICE function, adapting input and output
# for our code.
import spiceypy as spice

# The function cspice_kclear in MATLAB does not receive an input.

# The function spice.kclear in Python does not receive an input.

# The function cspice_kclear in MATLAB does not give an output:

# The function spice.kclear in Python does not give an output:

def mat2py_kclear():

    spice.kclear()
