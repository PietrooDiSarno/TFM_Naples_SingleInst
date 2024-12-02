# This code implements the functions that call the SPICE function, adapting input and output
# for our code.
import copy

import spiceypy as spice
import numpy as np

# The function cspice_ckgpav in MATLAB can receive:
# - inst: [1,1] = size(inst); int32 = class(inst)
# - sclkdp: [1,n] = size(sclkdp); double = class(sclkdp)
# - tol:[1,1] = size(tol); double = class(tol)
# - ref:[1,c1] = size(ref); char = class(ref)

# The function spice.ckgpav in Python receives:
# - inst: int
# - sclkdp: float --> a "for" cycle is required since MATLAB can handle over more values of "sclkdp", while Python cannot
# - tol: float OR int
# - ref: str

# The function cspice_ckgpav in MATLAB gives as output:
# - cmat: If [1,1] = size(sclkdp) then [3,3]   = size(cmat);double = class(cmat).
#         If [1,n] = size(sclkdp) then [3,3,n] = size(cmat); double = class(cmat)
# - av:[3,n] = size(av); double = class(av)
# - clkout:[1,n] = size(clkout); double = class(clkout)
# - found:[1,n] = size(found); logical = class(found)

# The function spice.ckgpav in Python gives as output:
# - tuple of (cmat,av,clkout,found). The type of the tuple is Tuple[ndarray, ndarray, float, bool]

def mat2py_ckgpav(inst,sclkdp,tol,ref):

    if np.size(sclkdp)==1:
        try:
          cmat,av,clkout=spice.ckgpav(inst,sclkdp,tol,ref)
          av=av.reshape(3,)
          found=True
        except Exception as e:
            if str(e) == 'Spice returns not found for function: ckgpav':
                cmat=np.zeros((3,3))
                av=np.zeros((3,))
                clkout=0
                found=False

    else:
      sclkdp = np.array(sclkdp).reshape(len(sclkdp), )
      n = len(sclkdp)
      cmat = np.array([])
      av = np.array([])
      clkout = np.array([])
      found = np.array([],dtype=bool)
      for i in range(n):
       try:
        Cmat,Av,Clkout=spice.ckgpav(inst,sclkdp[i],tol,ref)
        Av = Av.reshape(3, )
        if i == 0:
            av = copy.deepcopy(Av)
            cmat = copy.deepcopy(Cmat)
        else:
            cmat=np.stack((cmat,Cmat),axis=2)
            av=np.stack((av,Av),axis=1)
        clkout=np.append(clkout,Clkout)
        found=np.append(found,True)
       except Exception as e:
           if str(e) == 'Spice returns not found for function: ckgpav':
               if i == 0:
                   cmat = np.zeros((3, 3))
                   av = np.zeros((3, ))
               else:
                   cmat=np.stack((cmat,np.zeros((3, 3))),axis=2)
                   av=np.stack((av,np.zeros((3, ))),axis=1)
               clkout=np.append(clkout,0)
               found=np.append(found,False)

    return cmat,av,clkout,found
