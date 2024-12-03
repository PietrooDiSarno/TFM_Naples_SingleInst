# This code implements the functions that call the SPICE function, adapting input and output
# for our code.
import spiceypy as spice
import numpy as np

# The function cspice_getfov in MATLAB can receive:
# - instid: [1,1] = size(instid); int32 = class(instid)
# - room: [1,1] = size(room); int32 = class(room)

# The function spice.getfov in Python receives:
# - instid: int
# - room: int
# - shapelen: int (default value = 256)
# - framelen: int (default value = 256)

# The function cspice_getfov in MATLAB gives as output:
# - shape: [1,c1] = size(shape); char = class(shape)
# - frame: [1,c2] = size(frame); char = class(frame)
# - bsight: [3,1] = size(bsight); double = class(bsight)
# - bounds: [3,m] = [3,n] = size(bounds); double = class(bounds)

# The function spice.getfov in Python gives as output:
# - tuple of (shape,frame,bsight,n,bounds). The type of the tuple is Tuple[str, str, ndarray, int, ndarray]

def mat2py_getfov(instid,room):

    shape,frame,bsight,n,bounds_trans=spice.getfov(instid,room)

    bsight=bsight.reshape(3,)
    bounds = np.array([[row[i] for row in bounds_trans] for i in range(len(bounds_trans[0]))])

    return shape,frame,bsight,bounds