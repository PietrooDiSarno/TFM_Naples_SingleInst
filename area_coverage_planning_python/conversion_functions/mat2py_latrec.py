# This code implements the functions that call the SPICE function, adapting input and output
# for our code.
import numpy as np
import spiceypy as spice

# The function cspice_latrec in MATLAB can receive:
# - radius: [1,n] = size(radius); double = class(radius)
# - lon: [1,n] = size(lon); double = class(lon)
# - lat: [1,n] = size(lat); double = class(lat)

# The function spice.latrec in Python receives:
# - radius: float
# - lon: float
# - lat: float

# The function cspice_latrec in MATLAB gives as output:
# - rectan: [3,n] = size(rectan); double = class(rectan)

# The function spice.latrec in Python gives as output:
# - rectan: ndarray of shape (3,), "float" type

def mat2py_latrec(radius,lon,lat):

    if type(radius)==float:
        rectan=spice.latrec(radius,float(lon),float(lat))
        rectan=rectan.reshape(3,)

    else:
        radius=np.array(radius).reshape(np.size(radius),)
        lon=np.array(lon).reshape(np.size(lon),)
        lat=np.array(lat).reshape(np.size(lat),)
        rectan=np.zeros([3,np.size(radius)])
        for i in range(radius.size):
            rectan[:,i]=spice.latrec(float(radius[i]),float(lon[i]),float(lat[i]))

    return rectan
