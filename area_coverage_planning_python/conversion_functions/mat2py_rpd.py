# This code implements the functions that call the SPICE function, adapting input and output
# for our code.
import spiceypy as spice

# The function cspice_rpd in MATLAB does not receive an input.

# The function spice.rpd in Python does not receive an input.

# The function cspice_rpd in MATLAB gives as output:
# - rpd: DOUBLE = Scalar --> the number of radians per degree

# The function spice.rpd in Python gives as output:
# - rpd: float --> the number of radians per degree

def mat2py_rpd():
    rpd=spice.rpd()
    return rpd
