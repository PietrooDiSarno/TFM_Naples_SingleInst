# This code implements a function that receives the users' inputs, calls the SPICE function "spice.bodvrd"
# and returns the adapted outputs.
import spiceypy as spice
import numpy as np

# The function cspice_bodvrd in MATLAB can receive:
# - bodynm: [1,c1] = size(bodynm); char = class(bodynm) OR [1,1] = size(bodynm); cell = class(bodynm)
# - item: [1,c2] = size(item); char = class(item) OR [1,1] = size(item); cell = class(item)
# - maxn: [1,1] = size(maxn); int32 = class(maxn)

# The function spice.bodvrd receives,respectively:
# - bodynm: string
# - item: string
# - maxn: int

# When "bodynm" is given as a class "char" (1st case), it is expected to be a string in Python; when "bodynm" is
# given as a class "cell" in MATLAB (2nd case), in Python it is expected to be a list of one string. The same
# holds for "item".

# The function cspice_bodvrd in MATLAB gives as output:
# - values: [1,n] = size(values); double = class(values). (Trying this function in MATLAB, the real shape of the output
# is [n,1], so there is a reshape).
# The function spice.bodvrd in Python gives as output:
# - tuple of (dim, values). The type of the tuple is Tuple[int, ndarray]

# In order to make Python's output consistent with MATLAB's, in mat2py_bodvrd a while cycle is required:
# 'values' 's size in Python is equal to maxn, while in MATLAB it is equal to the number of values found (which
# corresponds to 'dim' in Python)


def mat2py_bodvrd(bodynm, item, maxn):

    if isinstance(bodynm,list): bodynm=bodynm[0]
    if isinstance(item,list): item=item[0]

    dim,values = spice.bodvrd(bodynm,item,maxn)
    while len(values)!=dim:
        values=np.delete(values,-1)
    values=np.reshape(values,[dim,])
    return values


