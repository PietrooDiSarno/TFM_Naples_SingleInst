from mosaic_algorithms.online_frontier_repair.frontierRepair import frontierRepair
from mosaic_algorithms.auxiliar_functions.plot.plotTour import plotTour
from mosaic_algorithms.paper.figure1.post_process_fig1 import post_process_fig1
from mosaic_algorithms.paper.figure1.input_data_fig1 import *  # Load mission info (kernels, SPICE ids, etc.)
import matplotlib.pyplot as plt

# Revision of grid discretization:
# Grid is going to be built in the camera frame, instead of the body-fixed
# frame. This is somewhat more complicated (it requires more calculations)
# but it could correct the spatial aberration that we presently see
# Re-implementation of sidewinder function with these new feature (and
# other code improvements)

# Online Frontier
A, fpList = frontierRepair(roistruct[0]['inittime'], stoptime, tcadence, inst, sc, target, roi, olapx, olapy, 3 * 1e-3)

# Plot tour
_ = plotTour(A, fpList, roistruct, sc, target)

# FOM post-process
post_process_fig1()

# Save figure [PDF]
figpath = '.'
plt.gcf().set_size_inches(9.2, 6)
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
filename = f"{figpath}/fig1_onlinefrontier.pdf"
plt.savefig(filename, dpi=1200, bbox_inches='tight')
