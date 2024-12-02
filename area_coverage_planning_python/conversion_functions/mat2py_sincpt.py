# This code implements a function that receives the users' inputs, calls the SPICE function "spice.sincpt"
# and returns the adapted outputs.

# The function cspice_sincpt in MATLAB can receive:
# - method: [1,c1] = size(method); char = class(method) OR [1,1] = size(method); cell = class(method)
# - target: [1,c2] = size(target); char = class(target) OR  [1,1] = size(target); cell = class(target)
# - et :  [1,1] = size(et); double = class(et)
# - fixref : [1,c3] = size(fixref); char = class(fixref) OR  [1,1] = size(fixref); cell = class(fixref)
# - abcorr: [1,c4] = size(abcorr); char = class(abcorr)  OR  [1,1] = size(abcorr); cell = class(abcorr)
# - obsrvr: [1,c5] = size(obsrvr); char = class(obsrvr) OR  [1,1] = size(obsrvr); cell = class(obsrvr)
# - dref: [1,c6] = size(dref); char = class(dref) OR [1,1] = size(dref); cell = class(dref)
# - dvec: [3,1] = size(dvec); double = class(dvec)

# The function spice.sincpt in Python receives:
# - method: string
# - target: string
# - et: float
# - fixref: string
# - abcorr: string
# - obsrvr: string
# - dref: string
# - dvec: numpy.ndarray of shape [3,1]

# The function cspice_sincpt in MATLAB gives as output:
# - spoint :  [3,1] = size(spoint); double = class(spoint)
# - trgepc:  [1,1] = size(trgepc); double = class(trgepc)
# - srfvec:[3,1] = size(srfvec); double = class(srfvec)
# - found:  [1,1] = size(found); logical = class(found)

# The function spice.sincpt in Python gives as output:
# - Tuple of [spoint,trgepc,srfvec,found] OR Tuple of [spoint,trgepc,srfvec] whose shape is
#   Tuple of [ndarray, float, ndarray, bool] or Tuple[ndarray, float, ndarray]

import spiceypy as spice
import numpy as np

def mat2py_sincpt(method,target,et,fixref,abcorr,obsrvr,dref,dvec):
    if isinstance(method,list): method=method[0]
    if isinstance(target,list): target=target[0]
    et=float(et)
    if isinstance(fixref,list): fixref=fixref[0]
    if isinstance(abcorr,list): abcorr=abcorr[0]
    if isinstance(obsrvr,list): obsrvr=obsrvr[0]
    if isinstance(dref,list): dref=dref[0]
    dvec=np.array(dvec).reshape(3,)
    try:
        spoint, trgepc, srfvec = spice.sincpt(method, target, et, fixref, abcorr, obsrvr, dref, dvec)
        spoint=np.array(spoint).reshape(3,)
        srfvec=np.array(srfvec).reshape(3,)
        found=True

    except  Exception as e:
         if str(e)=='Spice returns not found for function: sincpt':
           spoint=np.zeros([3,])
           trgepc=et
           srfvec=np.zeros([3,])
           found=False
    return spoint, trgepc, srfvec, found

