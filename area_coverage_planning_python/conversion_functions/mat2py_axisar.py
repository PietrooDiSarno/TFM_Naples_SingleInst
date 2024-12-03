# This code implements the functions that call the SPICE function, adapting input and output
# for our code.
import spiceypy as spice
import numpy as np

# The function cspice_axisar in MATLAB can receive:
# - axis: DOUBLE = Array[3,1]
# - angle: DOUBLE = Scalar

# The function spice.axisar in Python receives:
# - axis: ndarray(3,) or Iterable(float)
# - angle: float

# The function cspice_axisar in MATLAB gives as output:
# - r: DOUBLE=Array[3,3]

# The function spice.axisar in Python gives as output:
# - r: ndarray(3x3) of floats

def mat2py_axisar(axis, angle):
    axis=np.array(axis).reshape(3,)
    r=spice.axisar(axis,angle)

    return r
