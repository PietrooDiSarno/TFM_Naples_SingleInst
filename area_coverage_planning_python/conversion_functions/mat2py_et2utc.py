# This code implements the functions that call the SPICE function, adapting input and output
# for our code.
import spiceypy as spice
import numpy as np

# The function cspice_et2utc in MATLAB can receive:
# - et : [1,n] = size(et); double = class(et)
# - format: [1,c1] = size(format); char = class(format) OR  [1,1] = size(format); cell = class(format)
# - prec: [1,1] = size(prec); int32 = class(prec)

# The function spice.et2utc in Python receives:
# - et: float or iterable(float)
# - format: str
# - prec: int
# - lenout: int


# The function cspice_et2utc in MATLAB gives as output:
# - utcstr : [1,c2] = size(utcstr); char = class(utcstr)

# The function spice.et2utc in Python gives as output:
# - utcstr:ndarray or str

def mat2py_et2utc(et,format,prec):

    if isinstance(format,list):
        format=format[0]

    utcstr=spice.et2utc(et,format,prec)
    if np.size(et) > 1:
        utcstr=np.array(utcstr).reshape(len(et),)

    return utcstr
