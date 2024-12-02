import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from PIL import Image, ImageOps

def mapPlot(filename):

    # Future check: filename does not exist!

    root = tk.Tk()
    scrsz = root.winfo_screenwidth(), root.winfo_screenheight()
    root.withdraw()
    xtick = np.arange(-180, 181, 45)
    ytick = np.arange(-90, 91, 30)
    # x tick label
    xtickstr = []
    for x in xtick:
        if x < 0 and x > -180:
            xtickstr.append(f"{-x}ºW")
        elif x > 0 and x < 180:
            xtickstr.append(f"{x}ºE")
        else:
            xtickstr.append(f"{abs(x)}º")

    # y tick label
    ytickstr = []
    for y in ytick:
        if y < 0:
            ytickstr.append(f"{-y}ºS")
        elif y > 0:
            ytickstr.append(f"{y}ºN")
        else:
            ytickstr.append(f"{y}º")

    # figure
    plt.ion()
    W, L = scrsz
    fig, ax = plt.subplots(figsize=(W * 0.6 / 96, L * 0.5 / 96))
    #fig.subplots_adjust(left=0.4307, bottom=0.3144)
    map = Image.open(filename)
    #map = ImageOps.flip(map)
    ax.imshow(map, extent=(-180, 180, -90, 90), alpha=0.8, cmap='gray')
    ax.set_xticks(xtick)
    ax.set_yticks(ytick)
    ax.set_xticklabels(xtickstr)
    ax.set_yticklabels(ytickstr)
    ax.grid(True, color='w', linestyle=':', linewidth=0.5, alpha=1)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim([min(xtick), max(xtick)])
    ax.set_ylim([min(ytick), max(ytick)])
    ax.set_xlabel('Planetocentric longitude',fontsize = 20)
    ax.set_ylabel('Planetocentric latitude',fontsize = 20)
    ax.tick_params(axis = 'x', which = 'major', pad = 10)
    ax.tick_params(axis = 'y', which = 'major', pad = 10)
    ax.tick_params(labelsize=20)
    plt.pause(0.1)
    return fig, ax