# This code implements the functions that call the SPICE function, adapting input and output
# for our code.
import spiceypy as spice

# The function cspice_furnsh in MATLAB can receive:
# - file: [n,c1] = size(file); char = class(file) OR [1,n] = size(file); cell = class(file)

# The function spice.furnsh in Python receives:
# - file: str or iterable(str)

# The function cspice_furnsh in MATLAB gives no output.

# The function spice.furnsh in Python gives no output:

def mat2py_furnsh(file):

    if isinstance(file,str):
        spice.furnsh(file)

    else:
      for _, fl in enumerate(file):
        spice.furnsh(fl)
