# This code implements the functions that call the SPICE function, adapting input and output
# for our code.
import spiceypy as spice
import numpy as np

# The function cspice_tangpt in MATLAB can receive:
# - method:  [1,c1] = size(method); char = class(method) OR [1,1] = size(method); cell = class(method)
# - target: [1,c2] = size(target); char = class(target) OR [1,1] = size(target); cell = class(target)
# - et : [1,1] = size(et); double = class(et)
# - fixref: [1,c3] = size(fixref); char = class(fixref) OR [1,1] = size(fixref); cell = class(fixref)
# - abcorr: [1,c4] = size(abcorr); char = class(abcorr) OR [1,1] = size(abcorr); cell = class(abcorr)
# - corloc: [1,c5] = size(corloc); char = class(corloc) OR [1,1] = size(corloc); cell = class(corloc)
# - obsrvr: [1,c6] = size(obsrvr); char = class(obsrvr) OR [1,1] = size(obsrvr); cell = class(obsrvr)
# - dref: [1,c7] = size(dref); char = class(dref) OR [1,1] = size(dref); cell = class(dref)
# - dvec:[3,1] = size(dvec); double = class(dvec)

# The function spice.tangpt in Python receives:
# - method: str
# - target: str
# - et: float
# - fixref: str
# - abcorr: str
# - corloc: str
# - obsrvr: str
# - dref: str
# - dvec: ndarray or iterable(float)

# The function cspice_tangpt in MATLAB gives as output:
# - tanpt: [3,1] = size(tanpt); double = class(tanpt)
# - alt: [1,1] = size(alt); double = class(alt)
# - range: [1,1] = size(range); double = class(range)
# - srfpt: [3,1] = size(srfpt); double = class(srfpt)
# - trgepc: [1,1] = size(trgepc); double = class(trgepc)
# - srfvec: [3,1] = size(srfvec); double = class(srfvec)

# The function spice.tangpt in Python gives as output:
# - tuple of (tanpt,alt,range,srfpt,trgepc,srfvec). The type of the tuple is
# Tuple[ndarray,float,float,ndarray,float,ndarray]

def mat2py_tangpt(method,target,et,fixref,abcorr,corloc,obsrvr,dref,dvec):
    if isinstance(method,list): method=method[0]
    if isinstance(target,list): target=target[0]
    et=float(et)
    if isinstance(fixref,list): fixref=fixref[0]
    if isinstance(abcorr,list): abcorr=abcorr[0]
    if isinstance(obsrvr,list): obsrvr=obsrvr[0]
    if isinstance(dref,list): dref=dref[0]
    dvec=np.array(dvec).reshape(len(dvec),)

    tanpt,alt,range,srfpt,trgepc,srfvec=spice.tangpt(method,target,et,fixref,abcorr,corloc,obsrvr,dref,dvec)
    tanpt=np.array(tanpt).reshape(len(tanpt),)
    srfpt=np.array(srfpt).reshape(len(srfpt),)
    srfvec=np.array(srfvec).reshape(len(srfvec),)

    return tanpt,alt,range,srfpt,trgepc,srfvec
