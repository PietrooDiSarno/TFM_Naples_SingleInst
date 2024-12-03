# This code implements the functions that call the SPICE function, adapting input and output
# for our code.
import spiceypy as spice
import numpy as np

# The function cspice_sce2c in MATLAB can receive:
# - sc: [1,1] = size(sc); int32 = class(sc)
# - et: [1,n] = size(et); double = class(et)

# Since MATLAB can receive "et" also as a vector containing more values, while in Python this is not possible, a "for"
# cycle is required.

# The function spice.sce2c in Python receives:
# - sc: int
# - et: float

# The function cspice_sce2c in MATLAB gives as output:
# - sclkdp :[1,n] = size(sclkdp); double = class(sclkdp)

# The function spice.sce2c in Python gives as output:
# - sclkdp: float

def mat2py_sce2c(sc,et):

    if type(et)==float:
        sclkdp=spice.sce2c(sc,et)

    else:
       et = np.array(et).reshape(len(et), )
       sclkdp=[]
       for i in range(len(et)):
         sclkdp.append(spice.sce2c(sc,float(et[i])))
       sclkdp=np.array(sclkdp).reshape(len(sclkdp),)

    return sclkdp
