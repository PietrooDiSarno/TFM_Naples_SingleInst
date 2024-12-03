# This code implements the functions that call the SPICE function, adapting input and output
# for our code.
import numpy as np
import spiceypy as spice

# The function cspice_spkpos in MATLAB can receive:
# - normal: DOUBLE = Array[3]
# - point: DOUBLE = Array[3]

# The function spice.spkpos in Python receives:
# - normal: np.ndarray or iterable(floats)
# - point: np.ndarray or iterable(floats)

# The function cspice_spkpos in MATLAB gives as output:
# - plane: STRUCT = CSPICE_PLANE

# The function spice.spkpos in Python gives as output:
# - plane: plane

def mat2py_nvp2pl(normal,point):

    normal=np.array(normal,dtype='float').reshape(3,)
    point=np.array(point,dtype='float').reshape(3,)

    return spice.nvp2pl(normal,point)

