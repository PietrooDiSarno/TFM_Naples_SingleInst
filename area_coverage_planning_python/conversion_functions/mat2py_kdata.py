# This code implements the functions that call the SPICE functions, adapting input and output
# for our code.
import spiceypy as spice
import numpy as np

# The function cspice_kdata in MATLAB can receive:
# - which: [1,1] = size(which); int32 = class(which)
# - kind:  [1,c1] = size(kind); char = class(kind) OR [1,1] = size(kind); cell = class(kind)

# The function spice.kdata in Python receives:
# - which: int
# - kind: str
# - fillen: int; if not inserted, it is equal to the default value 256
# - typlen: int; if not inserted, it is equal to the default value 256
# - srclen: int; if not inserted, it is equal to the default value 256

# The function cspice_kdata in MATLAB gives as output:
# - file: [1,c2] = size(file); char = class(file)
# - filtyp: [1,c3] = size(filtyp); char = class(filtyp)
# - srcfil: [1,c4] = size(srcfil); char = class(srcfil)
# - handle: [1,1] = size(handle); int32 = class(handle)
# - found: [1,1] = size(found); logical = class(found)

# The function spice.kdata in Python gives as output:
# - tuple of (file,filtyp,srcfil,handle,found). The type of the tuple is Tuple[str, str, str, int, bool]

def mat2py_kdata(which,kind):

    if isinstance(kind,list):
        kind=kind[0]

    try:
      file,filtyp,srcfil,handle=spice.kdata(which,kind)
      found=True

    except Exception as e:
      if str(e)=='Spice returns not found for function: kdata':
          file=''
          filtyp=''
          srcfil=''
          handle=0
          found=False

    return file,filtyp,srcfil,handle,found
