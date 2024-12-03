# This code implements a function that receives the users' inputs, calls the SPICE function "pxform" and gives the
# adapted outputs
import spiceypy as spice
import numpy as np


# The function cspice_pxform in MATLAB can receive:
# - from: STRING = Scalar
# - to: STRING = Scalar
# - et: DOUBLE = Scalar  OR  DOUBLE = Array[N]

# The function spice.srfrec receives,respectively:
# - from: string
# - to: string
# - et: float

# The function cspice_pxform in MATLAB gives as output:
# - rotate: DOUBLE = Array[3,3]   OR   DOUBLE = Array[3,3,N]
# The function spice.pxform in Python gives as output:
# - rotate: ndarray

def mat2py_pxform(frm,to,et):

    if np.size(et)==1:
        rotate=spice.pxform(frm,to,et)

    else:
        et = np.array(et).reshape(len(et), )
        rotate = np.ndarray([3,3,len(et)],dtype=float)

        for i in range(et.size):
            rotate[:,:, i]= spice.pxform(frm,to,float(et[i]))

    return rotate


