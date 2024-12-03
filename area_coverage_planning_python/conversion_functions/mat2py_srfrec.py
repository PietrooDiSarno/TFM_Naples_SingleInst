# This code implements a function that receives the users' inputs, calls the SPICE function "spice.srfrec"
# and returns the adapted outputs.

import spiceypy as spice
import numpy as np

# The function cspice_srfrec in MATLAB can receive:
# - body:  [1,1] = size(body); int32 = class(body)
# - lon:  [1,n] = size(lon); double = class(lon)
# - lat: [1,n] = size(lat); double = class(lat)

# The function spice.srfrec in Python receives,respectively:
# - body: int
# - lon: float
# - lat: float

# When "lon" and "lat" are given as arrays of shape (1,n) in MATLAB, they are expected to be 1-D arrays in Python.
# Since Python cannot handle more than one value at a time both for "lon" and "lat", a "for" cycle is required.

# The function cspice_srfrec in MATLAB gives as output:
# - rectan: [3,n] = size(rectan); double = class(rectan)

# The function spice.srfrec in Python gives as output:
# - rectan: ndarray


def mat2py_srfrec(body,lon,lat):
    body=int(body)
    if np.size(lon)==1:
        lon=np.array(lon)
        lat=np.array(lat)
        rectan=spice.srfrec(body,float(lon),float(lat)).reshape(3,)

    else:
        rectan = np.zeros(shape=(3, np.size(lon)))
        lon=np.array(lon).reshape(np.size(lon),)
        lat=np.array(lat).reshape(np.size(lat),)
        for i in range(lon.size):
            rectan[:, i]= spice.srfrec(body,float(lon[i]),float(lat[i]))
    return rectan