# This code implements the functions that call the SPICE function, adapting input and output
# for our code.
import spiceypy as spice
import numpy as np

# The function cspice_fovray in MATLAB can receive:
# - inst:   [1,c1] = size(inst); char = class(inst) OR  [1,1] = size(inst); cell = class(inst)
# - raydir: [3,1] = size(raydir), double = class(raydir)
# - rframe: [1,c2] = size(rframe), char = class(rframe) OR  [1,1] = size(rframe); cell = class(rframe)
# - abcorr:  [1,c3] = size(abcorr), char = class(abcorr) OR [1,1] = size(abcorr); cell = class(abcorr)
# - obsrvr:  [1,c4] = size(obsrvr), char = class(obsrvr) OR  [1,1] = size(obsrvr); cell = class(obsrvr)
# - et: [1,n] = size(et), double = class(et)

# The function spice.fovray in Python receives:
# - inst: str
# - raydir: ndarray or Iterable(float)
# - rframe: str
# - abcorr: str
# - obsrvr: str
# - et:float

# The function cspice_bodn2c in MATLAB gives as output:
# - visibl:  [1,n] = size(visibl), logical = class(visibl)

# The function spice.bodn2c in Python gives as output:
# - visibl: bool

# Since MATLAB function can handle more values in "et", while Python cannot, a "for" cycle is required in order to make
# consistent inputs and outputs.

def mat2py_fovray(inst,raydir,rframe,abcorr,obsrvr,et):

    if isinstance(inst,list): bodynm=inst[0]
    if isinstance(rframe,list): rframe=rframe[0]
    if isinstance(abcorr,list): abcorr=abcorr[0]
    if isinstance(obsrvr,list): obsrvr=obsrvr[0]

    visibl=[]
    if np.size(et)==1:
        visibl.append(spice.fovray(inst, raydir, rframe, abcorr, obsrvr, et))
        return visibl
    for _,t in enumerate(et):
      visibl.append(spice.fovray(inst,raydir,rframe,abcorr,obsrvr,t))

    visibl = np.array(visibl)
    return visibl