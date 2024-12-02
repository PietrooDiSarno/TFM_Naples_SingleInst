# This code implements a function that receives the users' inputs, calls the SPICE function "spice.reclat"
# and returns the adapted outputs.
import numpy as np
import spiceypy as spice

# The function cspice_reclat in MATLAB can receive:
# - rectan: DOUBLE = Array[3]   or   DOUBLE = Array[3,N]

# The function spice.reclat in Python receives:
# - rectan: npy.ndarray or Iterable[float]

# The function cspice_reclat in MATLAB gives as output:
# - radius:  DOUBLE = Scalar   or   DOUBLE = Array[N]
# - lon: DOUBLE = Scalar   or   DOUBLE = Array[N]
# - lat: DOUBLE = Scalar   or   DOUBLE = Array[N]

# The function spice.reclat in Python gives as output:
# - tuple of (radius,lon,lat). The type of the tuple is Tuple[float,float,float]

def mat2py_reclat(rectan):
    rectan = np.array(rectan)

    if (rectan.size)==3:
        rectan=rectan.reshape(3,)
        result=spice.reclat(rectan)
        radius=result[0]
        lon=result[1]
        lat=result[2]
        return radius,lon,lat

    else:

       radius=np.zeros(rectan.shape[1])
       lon=np.zeros(rectan.shape[1])
       lat=np.zeros(rectan.shape[1])
       for i in range(rectan.shape[1]):
           rectan_array=np.array(rectan[:,i])
           result=spice.reclat(rectan_array)
           radius[i]=result[0]
           lon[i]=result[1]
           lat[i]=result[2]

    return radius,lon,lat


