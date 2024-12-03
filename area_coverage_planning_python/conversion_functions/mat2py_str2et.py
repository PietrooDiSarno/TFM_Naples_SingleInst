# This code implements the functions that call the SPICE function, adapting input and output
# for our code.

import spiceypy as spice
import numpy as np

# The function cspice_str2et in MATLAB can receive:
# - timstr: [n,c1] = size(timstr); char = class(timstr) OR [1,n] = size(timstr); cell = class(timstr)

# The function cspice_str2et in MATLAB can receive 2 types of inputs for "timstr": [n,c1]=size(timstr); char = class(timstr)
# OR  [1,n] =size(timstr); cell = class(timstr). We expect to receive a numpy.array of strings (shape : [n,1]) or a list
# of n strings.

# The function spice.str2et in Python receives:
# - timstr: str or iterable(str)

# The function cspice_str2et in MATLAB gives as output:
# - et: [1,n] = size(et); double = class(et)

# The function spice.str2et in Python gives as output:
# - et: float or ndarray(floats)

def mat2py_str2et(timstr):
    if isinstance(timstr,str):
        timstr=[timstr]
    if isinstance(timstr,np.ndarray):
        timstr=timstr.tolist()

    et = spice.str2et(timstr)
    et= np.array(et).reshape(np.size(et),)
    if len(et)==1:
        et=float(et[0])
    return et
