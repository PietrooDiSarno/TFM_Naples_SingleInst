import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

def post_process_fig3(roistruct,mosaic,index):

    # Zoom in
    plt.xlim([-135, -85])
    plt.ylim([-10, 40])
    xtick = range(-135, -84, 15)
    ytick = range(-10, 41, 10)

    degree_symbol = '$^\\circ$'

    # x tick label
    xtickstr = []
    for i in xtick:
        if i < 0 and i > -180:
            xtickstr.append(f"{-i}{degree_symbol}W")
        elif i > 0 and i < 180:
            xtickstr.append(f"{i}{degree_symbol}E")
        else:
            xtickstr.append(f"{abs(i)}{degree_symbol}")

    # y tick label
    ytickstr = []
    for i in ytick:
        if i < 0:
            ytickstr.append(f"{-i}{degree_symbol}S")
        elif i > 0:
            ytickstr.append(f"{i}{degree_symbol}N")
        else:
            ytickstr.append(f"{i}{degree_symbol}")

    plt.gca().set_xticks(xtick)
    plt.gca().set_yticks(ytick)
    plt.gca().set_xticklabels(xtickstr)
    plt.gca().set_yticklabels(ytickstr)
    plt.gca().tick_params(which='both', direction='in', top=True, right=True)
    plt.grid(True, which='both', color='w', linestyle=':', linewidth=1, alpha=1)
    plt.pause(3)

    # Save figure [PDF]
    figpath = 'png_images'
    plt.gcf().set_size_inches(6.96,5.5)
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
    roiname = roistruct[0]['name'].lower().replace(' ', '')
    name = f'post_process_{roiname}'
    filename = f"{figpath}/{roiname}_{mosaic}_{index}.png"
    plt.savefig(filename, dpi=1200, format='png', bbox_inches='tight')

