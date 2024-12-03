# This code implements the functions that call the SPICE function, adapting input and output
# for our code.
import spiceypy as spice

# The function cspice_bodn2c in MATLAB can receive:
# - name: STRING = Scalar

# The function spice.bodn2c in Python receives:
# - name: str

# The function cspice_bodn2c in MATLAB gives as output:
# - code: LONG = Scalar
# - found: BOOLEAN = Scalar

# The function spice.bodn2c in Python gives as output:
# - Tuple of [code,found] whose types are [int,bool] OR code:int

def mat2py_bodn2c(name):

    try:
     code=spice.bodn2c(name)
     found=True
    except Exception as e:
     if str(e) == 'Spice returns not found for function: bodn2c':
         code=0
         found=False

    return code,found