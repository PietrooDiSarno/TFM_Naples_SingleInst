import copy
import os
import pickle


from mosaic_algorithms.online_frontier_repair.frontierRepair import frontierRepair
from mosaic_algorithms.auxiliar_functions.plot.plotTour import plotTour
from mosaic_algorithms.paper.video.input_data_juice import *
import matplotlib.pyplot as plt
import importlib
import matplotlib

matplotlib.use('Agg')

mosaic = 'onlinefrontier'
roiname = roistruct[0]['name'].lower().replace(" ", "")
name = f'post_process_{roiname}'
module_name = f"mosaic_algorithms.paper.video.{name}"

folder = 'pickles_video'
for index, patron in enumerate(os.listdir(folder)):
    with open(f'{folder}/{patron}','rb') as f:
        A, fpList = pickle.load(f)
        # Plot tour
        ax = plotTour(A, fpList, roistruct, sc, target)
        ax.set_title(roistruct[0]['name'], fontweight = 'bold', fontsize = 20)
        """
        ax.legend(loc='upper center', ncol=2)
        handles, labels = ax.get_legend_handles_labels()
        vals_label = [labels[0]]
        inds = [0]
        for ind_label,val_label in enumerate(labels):
            if not val_label in vals_label:
                vals_label.append(val_label)
                inds.append(ind_label)
        labels = copy.deepcopy(vals_label)
        handles = [handles[i] for i in inds]
        handles = handles[:-1]
        labels = labels[:-1]
        legend = ax.legend(handles=handles, labels=labels, loc='upper center', ncol=2)
        legend.get_frame().set_alpha(1)
        """

        module = importlib.import_module(module_name)
        # FOM post-process
        exec(open(f"{name}.py").read())
        post_process_fig3(roistruct, mosaic, index)

