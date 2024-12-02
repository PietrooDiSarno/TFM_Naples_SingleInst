# This code implements a function that receives the users' inputs, calls the SPICE function "spice.subpnt"
# and returns the adapted outputs.
import spiceypy as spice
import numpy as np

# The function cspice_subpnt in MATLAB can receive:
# - method: [1,c1] = size(method); char = class(method) OR [1,1] = size(method); cell = class(method)
# - target: [1,c2] = size(target); char = class(target) OR [1,1] = size(target); cell = class(target)
# - et: [1,n] = size(et); double = class(et)
# - fixref: [1,c3] = size(fixref); char = class(fixref) OR [1,1] = size(fixref); cell = class(fixref)
# - abcorr: [1,c4] = size(abcorr); char = class(abcorr) OR [1,1] = size(abcorr); cell = class(abcorr)
# - obsrvr : [1,c5] = size(obsrvr); char = class(obsrvr) OR [1,1] = size(obsrvr); cell = class(obsrvr)

# The function spice.bodvrd receives,respectively:
# - method: string
# - target: string
# - et:  float
# - fixref: string
# - abcorr: string
# - obsrvr : string

# When "method" is given as a class "char" (1st case), it is expected to be a string in Python; when it is
# given as a class "cell" in MATLAB (2nd case), in Python it is expected to be a list of one string. The same
# holds for "target","fixref","abcorr","obsrvr".

# The function cspice_subpnt in MATLAB gives as output:
# - spoint: [3,n] = size(spoint); double = class(spoint)
# - trgepc: [1,n] = size(trgepc); double = class(trgepc)
# - srfvec: [3,n] = size(spoint); double = class(srfvec)
# The function spice.subpnt in Python gives as output:
# - tuple of (spoint,trgepc,srfvec). The type of the tuple is Tuple (ndarray, float, ndarray)

# Since MATLAB can give outputs for more values of "et", while Python can only handle one value at a time, a "for" cycle
# is required in mat2py_subpnt when n > 1: the outputs are charged in np.ndarrays("spoint" and "srfvec) and nd.array(trgepc)
# in order to make Python's output consistent with MATLAB's.

def mat2py_subpnt(method,target,et,fixref,abcorr,obsrvr):

    if isinstance(method, list): method = method[0]
    if isinstance(target, list): target = target[0]
    if isinstance(fixref, list): fixref = fixref[0]
    if isinstance(abcorr, list): abcorr = abcorr[0]
    if isinstance(obsrvr, list): obsrvr = obsrvr[0]

    if np.size(et)==1:
       spoint,trgepc,srfvec=spice.subpnt(method,target,et,fixref,abcorr,obsrvr)
       spoint=spoint.reshape(3,)
       srfvec=srfvec.reshape(3,)

    else:
        et=np.array(et).reshape(len(et),)
        spoint = np.zeros(shape=(3,et.size))
        trgepc = np.zeros(shape=(et.size,))
        srfvec = np.zeros(shape=(3, et.size))
        for i in range(et.size):
            spoint[:,i],trgepc[i],srfvec[:,i]=spice.subpnt(method,target,float(et[i]),fixref,abcorr,obsrvr)

    return spoint,trgepc,srfvec