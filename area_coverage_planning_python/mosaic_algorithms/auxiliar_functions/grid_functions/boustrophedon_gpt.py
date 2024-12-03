def boustrophedon(grid, dir1, dir2):
    """
    Plan a boustrophedon coverage path over a grid discretization.

    Parameters:
    -----------
    grid : list of lists
        A 2D grid represented as a list of lists (rows), where each element can be:
        - None: representing an empty cell.
        - A list or tuple [x, y]: representing coordinates of a point.

        Example:
        grid = [
            [None, [1, 2], None],
            [[3, 4], None, [5, 6]],
            [None, [7, 8], None]
        ]

    dir1 : str
        The primary sweeping direction. Must be one of:
        - 'north' or 'south' for horizontal sweeps.
        - 'east' or 'west' for vertical sweeps.

    dir2 : str
        The secondary sweeping direction orthogonal to dir1. Must be one of:
        - 'east' or 'west' if dir1 is 'north' or 'south'.
        - 'north' or 'south' if dir1 is 'east' or 'west'.

    Returns:
    --------
    tour : list
        A list of [x, y] coordinates representing the coverage tour.

    Notes:
    ------
    - Implements a boustrophedon (back-and-forth) path planning algorithm over the provided grid.
    - The sweeping starts from a corner determined by dir1 and dir2 and alternates direction after each row/column.
    - The function replicates the logic of the MATLAB version, adjusting for Python's 0-based indexing and data structures.

    Conversion Details:
    -------------------
    - MATLAB's 1-based indexing is converted to Python's 0-based indexing.
    - MATLAB cell arrays are replaced with Python lists.
    - MATLAB functions such as 'isempty', 'cellfun', and 'nnz' are replaced with Python equivalents.
    - The special case for grids with only one element is handled explicitly.
    """

    # Validate sweeping directions
    if dir1 in ['north', 'south']:
        if dir2 not in ['east', 'west']:
            raise ValueError("Sweeping direction is not well defined. dir2 must be 'east' or 'west' when dir1 is 'north' or 'south'.")
    elif dir1 in ['east', 'west']:
        if dir2 not in ['north', 'south']:
            raise ValueError("Sweeping direction is not well defined. dir2 must be 'north' or 'south' when dir1 is 'east' or 'west'.")
    else:
        raise ValueError("Invalid primary direction. dir1 must be 'north', 'south', 'east', or 'west'.")

    # Pre-allocate variables
    tour = []

    # Handle special case: grid with only one cell
    if len(grid) == 1 and len(grid[0]) == 1:
        cell = grid[0][0]
        if cell is not None:
            tour = [cell]  # Return the single cell in a list
        else:
            tour = []
        return tour

    # Determine initial sweep direction based on dir2
    if dir2 in ['east', 'south']:
        sweep = True
    else:
        sweep = False

    # Calculate the total number of non-empty cells to pre-allocate the tour list
    num_cells = sum(1 for row in grid for cell in row if cell is not None)
    tour = [None] * num_cells  # Pre-allocate list of planned observations

    ii = 0  # Index for tour

    # Plan tour over the grid discretization
    if dir1 in ['north', 'south']:  # Horizontal sweep
        if sweep:
            bearing = True  # left -> right
        else:
            bearing = False  # right -> left

        num_rows = len(grid)
        num_cols = len(grid[0])

        for i in range(num_rows):
            # Determine the row index based on the primary direction
            if dir1 == 'south':
                irow = i  # Sweep from top to bottom (0 to num_rows - 1)
            else:  # dir1 == 'north'
                irow = num_rows - i - 1  # Sweep from bottom to top

            # Safety check for index bounds (not strictly necessary here)
            if irow < 0 or irow >= num_rows:
                continue

            for j in range(num_cols):
                # Determine the column index based on the current bearing
                if not bearing:
                    icol = num_cols - j - 1  # right -> left
                else:
                    icol = j  # left -> right

                # Safety check for index bounds
                if icol < 0 or icol >= num_cols:
                    continue

                cell = grid[irow][icol]
                if cell is not None:
                    x, y = cell[0], cell[1]
                    tour[ii] = [x, y]  # Save the coordinate in the tour
                    ii += 1
            # Switch coverage direction after each row sweep
            bearing = not bearing

    elif dir1 in ['east', 'west']:  # Vertical sweep
        if sweep:
            bearing = True  # top -> down
        else:
            bearing = False  # bottom -> top

        num_rows = len(grid)
        num_cols = len(grid[0])

        for i in range(num_cols):
            # Determine the column index based on the primary direction
            if dir1 == 'west':
                icol = num_cols - i - 1  # Sweep from right to left
            else:  # dir1 == 'east'
                icol = i  # Sweep from left to right

            # Safety check for index bounds
            if icol < 0 or icol >= num_cols:
                continue

            for j in range(num_rows):
                # Determine the row index based on the current bearing
                if bearing:
                    irow = j  # top -> down
                else:
                    irow = num_rows - j - 1  # bottom -> top

                # Safety check for index bounds
                if irow < 0 or irow >= num_rows:
                    continue

                cell = grid[irow][icol]
                if cell is not None:
                    x, y = cell[0], cell[1]
                    tour[ii] = [x, y]  # Save the coordinate in the tour
                    ii += 1
            # Switch coverage direction after each column sweep
            bearing = not bearing
    else:
        # This case should not occur due to earlier validation, but included for completeness
        raise ValueError("Invalid primary direction. dir1 must be 'north', 'south', 'east', or 'west'.")

    return tour
