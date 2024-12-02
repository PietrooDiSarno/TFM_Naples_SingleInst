# This code implements the functions that call the SPICE function, adapting input and output
# for our code.
import spiceypy as spice
import numpy as np
import spiceypy.utils.support_types as stypes
import copy
# The function cspice_spkcov in MATLAB can receive:
# - spkfnm: [n,c1] = size(spkfnm); char = class(spkfnm) OR [1,n] = size(spkfnm); cell = class(spkfnm)
# - idcode: [1,1] = size(idcode); int32 = class(idcode)
# - room: [1,1] = size(room); int32 = class(room)
# - cover_i: [2m,1] = size(cover_i), double = class(cover_i) OR [0,0] = size(cover_i), double = class(cover_i) (optional)

# The function spice.spkcov in Python receives:
# - spkfnm: str
# - idcode: int
# - cover: SpiceCell(optional)

# The function cspice_spkcov in MATLAB gives as output:
# - cover: [2p,1] = size(cover), double = class(cover) OR [0,1] = size(cover), double = class(cover)

# The function spice.spkcov in Python gives as output:
# - cover: SpiceCell

def mat2py_spkcov(spkfnm,idcode,room,cover_i=None):

  if cover_i is not None:
    cover_i_vector=copy.deepcopy(cover_i)
    cover_i=stypes.SPICEDOUBLE_CELL(100)
    m=int(len(cover_i_vector)/2)
    for j in range(m):
      spice.wninsd(cover_i_vector[2*j],cover_i_vector[2*j+1],cover_i)


  if isinstance(spkfnm,list) or isinstance(spkfnm,np.ndarray):
     spkfnm=str(spkfnm[0])

  cover_cell = spice.spkcov(spkfnm,idcode,cover_i)
  cover=[]
  n=spice.wncard(cover_cell)

  for interval in range(n):
    start,end=spice.wnfetd(cover_cell,interval)
    cover.append(start)
    cover.append(end)

  cover=np.array(cover).reshape(np.size(cover),)
  return cover