# This code implements the functions that call the SPICE function, adapting input and output
# for our code.
import spiceypy as spice
import numpy as np

# The function cspice_srfnrm in MATLAB can receive:
# - method :  [1,c1] = size(method); char = class(method) OR  [1,1] = size(method); cell = class(method)
# - target: [1,c1] = size(target); char = class(target) OR  [1,1] = size(target); cell = class(target)
# - et: [1,1] = size(et); double = class(et)
# - fixref: [1,c1] = size(target); char = class(target) OR  [1,1] = size(target); cell = class(target)
# - srfpts:  [3,n] = size(srfpts); double = class(srfpts)

# The function spice.srfnrm in Python receives:
# - method: str
# - target: str
# - et: float
# - fixref: str
# - srfpts: ndarray of shape [n,3]

# The function cspice_srfnrm in MATLAB gives as output:
# - normls: [3,n] = size(normls); double = class(normls)

# The function spice.srfnrm in Python gives as output:
# - normls: ndarray of shape [n,3]

def mat2py_srfnrm(method,target,et,fixref,srfpts):

    if isinstance(method, list): method = method[0]
    if isinstance(target, list): target = target[0]
    if isinstance(fixref,list): fixref=fixref[0]
    if len(srfpts) == 3:
        srfpts=np.array(srfpts).reshape(1,3)
        normls = spice.srfnrm(method, target, et, fixref, srfpts) #output is array (1,3)
    else:
        srfpts_trans =np.array( [[row[i] for row in srfpts] for i in range(len(srfpts[0]))])
        normls_trans = spice.srfnrm(method,target,et,fixref,srfpts_trans)
        normls = np.array( [[row[i] for row in normls_trans] for i in range(len(normls_trans[0]))])
    return normls
