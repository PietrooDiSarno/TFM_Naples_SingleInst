# This code implements the function that calls the SPICE function, adapting input and output
# for our code.
import spiceypy as spice

# The function cspice_twopi in MATLAB does not receive an input.

# The function spice.twopi in Python does not receive an input.

# The function cspice_twopi in MATLAB gives as output:
# - twopi: [1,1] = size(twopi); double = class(twopi) --> the value of 2*pi

# The function spice.twopi in Python gives as output:
# - twopi: float --> the value of pi

def mat2py_twopi():

    twopi=spice.twopi()
    return twopi
