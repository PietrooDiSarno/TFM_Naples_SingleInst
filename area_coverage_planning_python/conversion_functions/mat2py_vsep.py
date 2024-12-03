# This code implements a function that receives the users' inputs, calls the SPICE function "spice.vsep"
# and returns the adapted outputs.
import spiceypy as spice
import numpy as np


# The function cspice_vsep in MATLAB can receive:
# - v1:[3,n] = size(v1); double = class(v1)
# - v2:[3,n] = size(v2); double = class(v2)
# Since MATLAB can handle over more than one couple of vectors, while Python can handle over one couple of vectors
# per time, a "for" cycle is needed.

# The function spice.vsep in Python receives:
# - v1: ndarray
# - v2: ndarray

# The function cspice_vsep in MATLAB gives as output:
# - vsep:[1,n] = size(vsep); double = class(vsep)

# The function spice.vsep in Python gives as output:
# - vsep:float

def mat2py_vsep(v1, v2):
    """
    Computes the vector separation between v1 and v2 using SPICE's vsep function.

    Warning: spice.vsep expects contiguous arrays in memory without any striding (i.e., the data should be laid out
    sequentially in memory). When you extract slices from a NumPy array using operations like v1[:, i], the resulting
    arrays may not be contiguous, leading to the TypeError: strided arrays not supported.
    To resolve this issue, you need to ensure that the vectors you pass to spice.vsep are contiguous in memory and
    have the appropriate data type (float64). You can achieve this by using numpy.ascontiguousarray or numpy.copy()
    with the correct data type.

    Parameters:
    - v1: ndarray of shape (3,) or (3, N)
    - v2: ndarray of shape (3,) or (3, N)

    Returns:
    - vsep: ndarray of shape (3,) if single vectors are provided,
            or (3, N) for multiple vectors.
    """

    # Copy input arrays as float64
    v1 = np.asarray(v1, dtype=np.float64)
    v2 = np.asarray(v2, dtype=np.float64)

    # Check if inputs are single vectors
    if v1.ndim == 1 and v1.size == 3 and v2.ndim == 1 and v2.size == 3:
        # Ensure the vectors are 1D arrays with 3 elements
        vsep = spice.vsep(v1, v2)

    # Check if inputs are arrays of multiple vectors
    elif v1.ndim == 2 and v1.shape[0] == 3 and v2.ndim == 2 and v2.shape[0] == 3:
        num_vectors = v1.shape[1]

        if v2.shape[1] != num_vectors:
            raise ValueError("v1 and v2 must have the same number of vectors along axis 1.")

        vsep_list = []
        for i in range(num_vectors):
            # Extract the i-th vector from each input
            vec1 = v1[:, i]
            vec2 = v2[:, i]

            # Ensure the vectors are contiguous and of type float64
            vec1_contig = np.ascontiguousarray(vec1, dtype=np.float64)
            vec2_contig = np.ascontiguousarray(vec2, dtype=np.float64)

            # Compute the vector separation
            vsep_result = spice.vsep(vec1_contig, vec2_contig)
            vsep_list.append(vsep_result)

        # Convert the list of vsep results back to a NumPy array
        vsep = np.column_stack(vsep_list)

    else:
        raise ValueError("Input vectors must be either 1D arrays of length 3 or 2D arrays with shape (3, N).")

    return vsep
