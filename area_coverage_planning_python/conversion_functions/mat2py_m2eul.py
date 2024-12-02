# This code implements the functions that call the SPICE function, adapting input and output
# for our code.
import spiceypy as spice
import numpy as np

# The function cspice_m2eul in MATLAB can receive:
# - r: DOUBLE = Array[3,3]
# - axis3: LONG = Scalar
# - axis2: LONG = Scalar
# - axis1: LONG = Scalar

# The function spice.m2eul in Python receives:
# - r: ndarray or iterable of iterables(float)
# - axis3: int
# - axis2: int
# - axis1: int

# The function cspice_m2eul in MATLAB gives as output:
# - angle3:  DOUBLE = Scalar
# - angle2:  DOUBLE = Scalar
# - angle1:  DOUBLE = Scalar

# The function spice.bodn2c in Python gives as output:
# - tuple of (angle3,angle2,angle1). The type of the tuple is Tuple[float,float,float]

def mat2py_m2eul(r,axis3,axis2,axis1):

    r=np.array(r).reshape(3,3)
    axis3=int(axis3)
    axis2=int(axis2)
    axis1=int(axis1)

    angle3,angle2,angle1=spice.m2eul(r,axis3,axis2,axis1)
    angle3=float(angle3)
    angle2=float(angle2)
    angle1=float(angle1)

    return angle3,angle2,angle1