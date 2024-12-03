import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

def post_process_fig1():
    # Zoom in
    plt.xlim([-30, 60])
    plt.ylim([-15, 35])
    xtick = range(-30, 61, 15)
    ytick = range(-15, 36, 10)

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