# This code implements the functions that call the SPICE function, adapting input and output
# for our code.
import numpy as np
import spiceypy as spice

# The function cspice_spkpos in MATLAB can receive:
# - targ: STRING = Scalar
# - et: DOUBLE = Scalar   OR   DOUBLE = Array[N]
# - ref: STRING = Scalar
# - abcorr: STRING = Scalar
# - obs: STRING = Scalar

# The function spice.spkpos in Python receives:
# - targ: string
# - et: float or ndarray
# - ref: string
# - abcorr: string
# - obs:string

# The function cspice_spkpos in MATLAB gives as output:
# - ptarg: DOUBLE = Array[3]   or   DOUBLE = Array[3,N]
# - ltime: DOUBLE = Scalar   or   DOUBLE = Array[N]

# The function spice.spkpos in Python gives as output:
# - tuple of (ptarg,ltime). The type of the tuple is Tuple[ndarray, float] or Tuple[ndarray,ndarray]

def mat2py_spkpos(targ,et,ref,abcorr,obs):
    if np.size(et) == 1:
        et=float(et)
        ptarg, ltime = spice.spkpos(targ, et, ref, abcorr, obs)
        ptarg=np.array(ptarg,'float').reshape(3,)
        ltime=float(ltime)
    else:
        et=np.array(et,dtype='float').reshape((len(et),))
        ptarg_trans, ltime = spice.spkpos(targ, et, ref, abcorr, obs)
        ptarg = np.array([[row[i] for row in ptarg_trans] for i in range(len(ptarg_trans[0]))],dtype='float')
        ltime=np.array(ltime,dtype='float').reshape((len(et),))
    return ptarg,ltime



