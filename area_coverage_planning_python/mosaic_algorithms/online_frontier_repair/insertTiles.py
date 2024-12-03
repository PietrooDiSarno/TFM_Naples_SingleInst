import copy
import numpy as np
from conversion_functions import *


def insertTiles(*args):
    """
    This function includes new observation points in a planned tour
    observations

    Programmers:  Paula Betriu (UPC/ESEIAAT)
    Date:         06/2023

    Usage:        [tour, map] = insertTiles(map, newp, indp)

    Inputs:
      > map:       list of lists of grid points. In order to avoid
                   mapping boundaries, map is bounded by NaN rows and
                   columns (first and last)
      > newp:      list of lists of the new observation points to be
                   included in 'tour'
      > indp:      list of lists of the new observation points locations to
                   be included in 'map'

    Returns:
      > map:       updated list of lists of grid points

    """

    map = args[0]
    newp = args[1]
    indp = args[2]

    # Insert elements in map
    offcol = 0
    offrow = 0

    for i in range(len(newp)):
        aux_map = copy.deepcopy(map)  # save the map at its current state (for potential relocation purposes)
        indel = indp[i]

        # Update index position
        indel[0] += offrow
        indel[1] += offcol

        # If the index is in the boundaries (or even further) of the map, then
        # we will have to relocate the elements in map and create a bigger grid

        offrow0 = copy.deepcopy(offrow)
        offcol0 = copy.deepcopy(offcol)

        if indel[0] >= (len(map) - 1):  # last row or more
            nrows = 2 + indel[0] - len(map)  # number of additional rows in the grid (last rows)

            # Map relocation
            map =  aux_map + [[np.array([np.nan,np.nan]) for _ in range(len(map[0]))] for _ in range(nrows)]

        elif indel[0] <= 0:  # First row or less
            nrows = 1 - indel[0]  # Number of additional rows in the grid (first rows)
            offrow += nrows # rows offset

            # Map relocation
            map = [[np.array([np.nan,np.nan]) for _ in range(len(map[0]))] for _ in range(nrows)] + aux_map

        aux_map = copy.deepcopy(map)

        if indel[1] >= (len(map[0]) -1) :  # last column or more
            ncols = 2 + indel[1] - len(map[0])  # number of additional columns
            # in the grid (last columns)

            # Map relocation
            map = [row + [np.array([np.nan, np.nan])] * ncols for row in aux_map]

        elif indel[1] <= 0:  # First column or less
            ncols = 1 - indel[1]  # number of additional columns in the grid (first columns)
            offcol += ncols # columns offset

            # Map relocation
            map = [[np.array([np.nan, np.nan])] * ncols + row for row in aux_map]

        # Update the current element index
        indel[0] += (offrow - offrow0)
        indel[1] += (offcol - offcol0)

        # Include elements in map
        map[indel[0]][indel[1]] = newp[i]

    return map

# The following is a MATLAB code that must be translated in Python. Since it is a comment, it is not used so it is not
# translated for now.

# [Future work]: # # Insert elements in tour?
# switch heuristics
#     case 'nearest neighbour'
#
#         # Find the nearest neighbour in the tour and insert the new
#         # observation before or after this element
#         for i=1:length(newp)
#
#             # Find the nearest element in the tour
#             mindist = inf;
#             for j=1:length(tour)
#                 dist = norm(tour{j} - newp{i});
#                 if dist < mindist
#                     mindist = dist;
#                     idx     = j;
#                 end
#             end
#
#             # Evaluate before and after elements distance to decide
#             # whether the new point should be put before or after its
#             # nearest neighbour
#             aux_tour = tour;
#             if idx ~= 1 && idx < length(tour)
#                 dist1 = norm(tour{idx - 1} - newp{i});
#                 dist2 = norm(tour{idx + 1} - newp{i});
#                 if dist1 < dist2
#                     tour(idx)     = newp(i);
#                     tour(idx+1:end+1) = aux_tour(idx:length(tour));
#                 else
#                     tour(idx+1)     = newp(i);
#                     tour(idx+2:length(tour)+1) = aux_tour(idx+1: ...
#                         length(tour));
#                 end
#
#             # Particular cases: when the nearest neighbour is the first
#             # or last element in the tour
#             elseif idx == 1
#                 dist = norm(tour{2} - newp{i});
#                 if dist < mindist
#                     tour(2) = newp(i);
#                     tour(3:length(tour)+1) = aux_tour(2:length(tour));
#                 else
#                     tour(1) = newp(i);
#                     tour(2:length(tour)+1) = aux_tour(:);
#                 end
#             else
#                 dist = norm(tour{idx - 1} - newp{i});
#                 if dist < mindist
#                     tour(idx)     = newp(i);
#                     tour(idx + 1) = aux_tour(idx);
#                 else
#                     tour(idx + 1) = newp(i);
#                 end
#             end
#         end
#     case 'manhattan'
#
#         # Analyze the tour length when inserting the point at each
#         # position in the tour
#         for i=1:length(newp)
#             mindelta = inf;
#             for j=1:length(tour)
#                 if j ~= length(tour)
#                     %oldd = norm(tour{j+1} - tour{j});
#                     oldd = manhattan(tour{j+1}, tour{j});
#                     %newd = norm(tour{j+1} - newp{i}) + ...
#                     %    norm(tour{j} - newp{i});
#                     newd = manhattan(tour{j+1}, newp{i}) + ...
#                         manhattan(tour{j}, newp{i});
#                     delta = newd - oldd; % increase of tour length when
#                     % adding the new insertion point
#                 else
#                     %delta = norm(tour{j} - newp{i});
#                     delta = manhattan(tour{j}, newp{i});
#                 end
#                 if delta < mindelta
#                     mindelta = delta;
#                     idx = j;
#                 end
#             end
#
#             # Insert element in the position that minimizes the tour
#             # deviation length
#             if idx ~= length(tour)
#                 aux_tour = tour;
#                 tour(idx+1) = newp(i);
#                 tour(idx+2:end+1) = aux_tour(idx+1:end);
#             else
#                 tour(end+1) = newp(i);
#             end
#
#         end
#     case 'simulated annealing'
#
# end