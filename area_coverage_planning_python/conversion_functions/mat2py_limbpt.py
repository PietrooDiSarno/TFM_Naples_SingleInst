# This code implements a function that receives the users' inputs, calls the SPICE function "spice.limbpt"
# and returns the adapted outputs.
import spiceypy as spice
import numpy as np

# The function cspice_limbpt in MATLAB can receive:
# - method: [1,c1] = size(method); char = class(method) OR [1,1] = size(method); cell = class(method)
# - target: [1,c2] = size(target); char = class(target) OR [1,1] = size(target); cell = class(target)
# - et:  [1,1] = size(et); double = class(et)
# - fixref: [1,c3] = size(fixref); char = class(fixref) OR [1,1] = size(fixref); cell = class(fixref)
# - abcorr: [1,c4] = size(abcorr); char = class(abcorr) OR [1,1] = size(abcorr); cell = class(abcorr)
# - corloc: [1,c5] = size(corloc); char = class(corloc) OR [1,1] = size(corloc); cell = class(corloc)
# - obsrvr: [1,c6] = size(obsrvr); char = class(obsrvr) OR [1,1] = size(obsrvr); cell = class(obsrvr)
# - refvec: [3,1] = size(refvec); double = class(refvec)
# - rolstp: [1,1] = size(rolstp); double = class(rolstp)
# - ncuts: [1,1] = size(ncuts); int32 = class(ncuts)
# - schstp: [1,1] = size(schstp); double = class(schstp)
# - soltol: [1,1] = size(soltol); double = class(soltol)
# - maxn: [1,1] = size(maxn); int32 = class(maxn)

# If "method" is given as a string (1st case) in MATLAB, it is expected to be a string in Python; if it is given as
# a cell (2nd case) whose only element is the string, it is expected to be a list containing only the string in Python;
# the same for "target", "fixref","abcorr","corloc","obsrvr".


# The function spice.limbpt in Python receives,respectively:
# - method: string
# - target: string
# - et: float
# - fixref: string
# - abcorr: string
# - corloc: string
# - obsrvr: string
# - refvec: ndarray OR iterable of floats
# - rolstp: float
# - ncuts: int
# - schstp: float
# - soltol: float
# - maxn: int

# The function cspice_limbpt in MATLAB gives as output:
# - npts: [1,n] = size(npts); int32 = class(npts)
# - points: [3,m] = size(points); double = class(points)
# - epochs: [1,m] = size(epochs); double = class(epochs)
# - tangts: [3,m] = size(tangts); double = class(tangts)

# The function spice.limbpt in Python gives as output:
# - tuple of (npts,points,epochs,tangts). The type of the tuple is Tuple[ndarray, ndarray,ndarray,ndarray]

def mat2py_limbpt(method,target,et,fixref,abcorr,corloc,obsrvr,refvec,rolstp,ncuts,schstp,soltol,maxn):
    if isinstance(method, list): method = method[0]
    if isinstance(target, list): target = target[0]
    if isinstance(fixref, list): fixref = fixref[0]
    if isinstance(abcorr, list): abcorr = abcorr[0]
    if isinstance(obsrvr, list): obsrvr = obsrvr[0]
    if isinstance(corloc,list): corloc = corloc[0]

    refvec=np.array(refvec).reshape(3,)

    npts,points_trans,epochs,tangts_trans = spice.limbpt(method,target,et,fixref,abcorr,corloc,obsrvr,refvec,rolstp,ncuts,schstp,soltol,maxn)
    nf=len(npts)

    points = np.array([[row[i] for row in points_trans] for i in range(len(points_trans[0]))])
    tangts = np.array([[row[i] for row in tangts_trans] for i in range(len(tangts_trans[0]))])

    if nf!=maxn:
        npts=np.concatenate((npts,np.zeros(maxn-nf)))
        points=np.hstack((points,np.zeros([3,maxn-nf])))
        epochs=np.concatenate((epochs,np.zeros(maxn-nf)))
        tangts = np.hstack((tangts, np.zeros([3,maxn-nf])))

    return npts,points,epochs,tangts

