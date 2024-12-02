# This code implements the functions that call the SPICE function, adapting input and output
# for our code.
import spiceypy as spice

# The function cspice_dpr in MATLAB does not receive an input.

# The function spice.dpr in Python does not receive an input.

# The function cspice_dpr in MATLAB gives as output:
# - dpr: [1,1] = size(dpr); double = class(dpr) --> the number of degrees per radian

# The function spice.dpr in Python gives as output:
# - dpr: float --> the number of degrees per radian

def mat2py_dpr():
    dpr=spice.dpr()
    return dpr

