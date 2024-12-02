# This code implements the function that call the SPICE function, adapting input and output
# for our code.
import spiceypy as spice
import numpy as np

# The function cspice_cnmfrm in MATLAB receives "cname" as input.
# It can receive 2 types of inputs for "cname":  [n,c1]=size(cname); char = class(cname)
# OR  [1,n] =size(cname); cell = class(cname). In Python we expect to receive a numpy.array of strings OR a list of
# strings; if n = 1, in the first case we expect to receive a string, in the second case we expect to receive
# a list containing a single string.

# The function spice.cnmfrm in Python receiven "cname" as input:
# cname: string
# lenout: int (default value = 256 )

# The function cnmfrm in MATLAB gives the following outputs:
# - frcode:[1,n]=size(frcode); int32=class(frcode)
# - frname: [n,c2]=size(frname); Char=class(frname)
# - found: [1,n] = size(found); Logical=class(found)

# The function spice.cnmfrm gives the following outputs:
# - Tuple of [frcode,frname,found] OR [frcode,frname] whose type is Tuple[int, str, bool] OR Tuple[int, str]

# mat2py_cnmfrm returns:
# - frcode: numpy.array of int whose shape is [n,] (1-D)
# - frname: list of n strings
# - found: numpy.array of bool whose shape is [n,] (1-D)


def mat2py_cnmfrm(cname):
    if isinstance(cname,str):
        cname=[cname]
    frcode = []
    frname = []
    found = []

    #spice.cnmfrm takes only one name as input (differently from cspice_cnmfrm):
    #we iterate over the names contained in cname
    for _,name in enumerate(cname):

      try:
        frcode.append(spice.cnmfrm(name)[0])
        frname.append(spice.cnmfrm(name)[1])
        found.append(1)
      except Exception as e:
         print(f'Exception for {name}: {e}')
         frname.append(' ')
         found.append(0)
         frcode.append(0)

    if len(frname)==1:
        return frcode[0],frname[0], found[0]
    else:
        frcode=np.array(frcode,dtype='int32')
        found=np.array(found,dtype='bool')

        return frcode,frname,found





