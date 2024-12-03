# This code implements the functions that call the SPICE function, adapting input and output
# for our code.
import spiceypy as spice

# The function cspice_ktotal in MATLAB can receive:
# - kind :  [1,c1] = size(kind); char = class(kind) OR [1,1] = size(kind); cell = class(kind)

# The function spice.ktotal in Python receives:
# - kind: str

# The function cspice_ktotal in MATLAB gives as output:
# - count: [1,1] = size(count); int32 = class(count)

# The function spice.ktotal in Python gives as output:
# - count: int

def mat2py_ktotal(kind):

    if isinstance(kind,list):
        kind=kind[0]

    count=spice.ktotal(kind)

    return count
