import matplotlib as plt
from datetime import datetime
import spiceypy as spice
import numpy as np
from FuturePackage import DataManager

class DataPlotter:
    def __init__(self, cParams):
        self.data = cParams

    def plot_plan(self, col=None):
        roiL = DataManager.getInstance().getRoiList()
        data = self.data
        if col is None:
            col = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
        color = dict(zip(roiL[:].name, col))
        H = 0.5  # just to plot tw, height of each roi tw
        cc = 0
        for tw in range(spice.wncard(self.searchSpace)):
            fig, ax = plt.subplots()
            _, s, e = self.getAllTW(tw)
            ftv = np.linspace(s, e, 6000)
            for i, roi in enumerate(self.roitwl):
                cc = cc + 1
                if cc > len(col) - 1: cc = 0
                q = []
                for ft in ftv:
                    q.append(self.evalQualityRoi(i, ft))
                ax.plot(ftv, q, label=roi.name, color=color[roi.name])
                for interval in range(spice.wncard(roi)):
                    start, end = spice.wnfetd(roi, interval)
                    if self.stol[i] >= start and self.stol[i] <= end:
                        tw = newTimeWindow(start, end)
                        plot_tw(ax, tw, 0 + i * H, 0 + H * (i + 1), color[self.roikeyl[i]])  # plot feasible tw
                        obstw = newTimeWindow(self.stol[i], self.stol[i] + self.obsLength[i])
                        plot_tw(ax, obstw, 0 + i * H, 0 + H * (i + 1), 'w')  # plot feasible tw
        fig_numbers = plt.get_fignums()
        figs_with_nodata = [plt.figure(num) for num in fig_numbers if not plt.figure(num).axes[0].has_data()]
        if len(figs_with_nodata) != 0:
            for fig in figs_with_nodata:
                plt.close(fig)
        plt.show()

    def plot_gantt(self, col=None):
        roi = self.roikeyl
        stol = self.stol
        # obslength = self.obslength
        sorted_index = sorted(range(len(stol)), key=lambda i: stol[i])
        stol = [stol[i] for i in sorted_index]
        roi = [roi[i] for i in sorted_index]
        # obslength = [obslength[i] for i in sorted_index]
        lo = np.linspace(0, 1, len(self.stol) + 1)
        f, axs = plt.subplots(1, len(self.stol), sharey=False, facecolor='w')
        for i in range(len(axs)):
            intbeg = stol[i]
            intend = intbeg + 800  # +obslength
            axs[i].add_patch(
                patches.Rectangle((intbeg, lo[-2 - i]), width=intend - intbeg, height=lo[-1 - i] - lo[-2 - i], lw=1,
                                  color='b', fill=True))
            axs[i].set_xlim(intbeg, intend)

            if i == 0:
                axs[i].spines['right'].set_visible(False)
                axs[i].yaxis.set_ticks([])

            elif i == len(axs) - 1:
                axs[i].spines['left'].set_visible(False)
            else:
                axs[i].spines['right'].set_visible(False)
                axs[i].spines['left'].set_visible(False)
            axs[i].set_xticks([intbeg, intend])
            str_beg = spice.et2utc(intbeg, 'C', 0)
            str_end = spice.et2utc(intend, 'C', 0)

            datetime_beg = datetime.strptime(str_beg, '%Y %b %d %H:%M:%S')
            datetime_end = datetime.strptime(str_end, '%Y %b %d %H:%M:%S')
            str_beg = datetime_beg.strftime('%d/%m/%Y')
            str_end = datetime_end.strftime('%d/%m/%Y')

            time_beg = datetime_beg.strftime('%H:%M:%S')
            time_end = datetime_end.strftime('%H:%M:%S')

            y_height = (lo[-1 - i] + lo[-2 - i]) / 2
            axs[i].set_yticks([y_height])
            roi_name = [roi[i]]
            axs[i].set_yticklabels(roi_name)
            axs[i].set_xticklabels([time_beg, time_end], ha='right', va='center')
            axs[i].tick_params(axis='x', rotation=45, pad=20)

            if str_beg == str_end:
                axs[i].annotate(str_beg, ((intbeg + intend) / 2, (lo[-1 - i] + lo[-2 - i]) / 2), color='white',
                                weight='bold', fontsize=10, ha='center', va='center')
            else:
                axs[i].annotate(str_beg + '-' + str_end, ((intbeg + intend) / 2, (lo[-1 - i] + lo[-2 - i]) / 2),
                                color='white', weight='bold', fontsize=10, ha='center', va='center')
        roi_names = []
        y_heights = []
        # for i in range(len(axs)):
        # height = (lo[-1-i]+lo[-2-i])/2
        # len()

        plt.title('Mission Gantt Chart')
        plt.show()

        # plt.tight_layout()

        # def generateObsDataBase(self, tw, step, roi, obs, ifov, npix, imageRate, fs):
        # tv = getSteps_tw(tw, step)
        # nimages = [0] * len(tv)
        # timeobs = [0] * len(tv)

        # for i, t in enumerate(tv):
        #    nimages[i], timeobs[i] = self.aproxObsTime(t, aroi1, 'GLL', ifov, npix, imageRate, fs)

        # return [tv, nimages, timeobs]

    def plotTWRoiParams(self, roi, index, detailed=False):
        if index == 2:
            datalabel = 'resolution [km/px]'
        else:
            datalabel = 'duration [s]'
        time = self.obsDataBase[self.roikeyl[roi]][0]
        data = self.obsDataBase[self.roikeyl[roi]][index]
        f, axs = plt.subplots(1, spice.wncard(self.roitwl[roi]), sharey=False, facecolor='w')
        for i in range(len(time)):
            t = time[i][:]
            d = data[i][:]
            if detailed:
                if index == 2:
                    mymin, col = self.findMin(d)
                    for j, elem in enumerate(d):
                        if elem <= d[-1]:
                            break
                else:
                    for j, elem in enumerate(d):
                        if elem > max(d) / 2:
                            break
                t = time[i][j:]
                d = data[i][j:]
            axs[i].plot(t, d)
            axs[i].set_xlim(t[0], t[-1])

            str_beg = spice.et2utc(t[0], 'C', 0)
            str_end = spice.et2utc(t[-1], 'C', 0)

            datetime_beg = datetime.strptime(str_beg, '%Y %b %d %H:%M:%S')
            datetime_end = datetime.strptime(str_end, '%Y %b %d %H:%M:%S')

            str_beg = datetime_beg.strftime('%d %b %Y')
            str_end = datetime_end.strftime('%d %b %Y')

            if str_beg == str_end:
                axs[i].set_title(str_beg)
            else:
                axs[i].set_title(str_beg + '-\n' + str_end)

            labelFormat = '%H:%M:%S'
            tlabel, labels = self.writeTWlabels(t, labelFormat, nlabels=5)

            axs[i].set_xticks(tlabel)
            axs[i].tick_params(axis='x', rotation=45, pad=20)
            axs[i].set_xticklabels(labels, ha='right', va='center')
        f.suptitle('Observation ' + datalabel + ' for: ' + self.roikeyl[roi])

    def writeTWlabels(self, time, format, nlabels=5):
        t = np.linspace(time[0], time[-1], nlabels)
        aux = []
        for h in range(len(t)):
            j = datetime.strptime(spice.et2utc(t[h], 'C', 0), '%Y %b %d %H:%M:%S')
            aux.append(j.strftime(format))
        return t, aux

    def findMin(self, data):
        mymin = 9999.
        col = -99999
        for k, elem in enumerate(data):
            if mymin > elem:
                mymin = elem
                col = k
        return mymin, col

    def plotMaxRes(self, roi, start=None, end=None):  # returns the time (seconds) needed to perform observation number i
        globmin = 99999.
        res = []
        time = self.obsDataBase[self.roikeyl[roi]][0]

        # card = spice.wncard(self.roitwl[roi])
        # start, _ = spice.wnfetd(self.roitwl[roi],0)
        # _,end = spice.wnfetd(self.roitwl[roi], card-1)
        if start is None and end is None:
            instRes = self.obsDataBase[self.roikeyl[roi]][2]
            for i, vect in enumerate(instRes):
                mymin, j = self.findMin(vect)
                if globmin > mymin:
                    globmin = mymin
                    col = j
                    row = i
            start = time[row][col] - 24 * 3600 * 5
            end = time[row][col] + 24 * 3600 * 5
            print(globmin)
            print(spice.et2utc(time[row][col], 'C', 0))
        t = np.linspace(start, end, 6000)
        startY = datetime.strptime(spice.et2utc(t[0], 'C', 0), '%Y %b %d %H:%M:%S')
        startY = startY.strftime('%Y')
        endY = datetime.strptime(spice.et2utc(t[-1], 'C', 0), '%Y %b %d %H:%M:%S')
        endY = endY.strftime('%Y')

        if startY != endY:
            titleY = startY + '-' + endY
        else:
            titleY = startY
        for i, et in enumerate(t):
            res.append(self.evalQualityRoi(roi, et))
            if np.isnan(res[i]):
                res[i] = 300.
        labelFormat = "%d %b"
        tplot, labels = self.writeTWlabels(t, labelFormat, nlabels=10)

        plt.figure()
        plt.plot(t, res)
        plt.xticks(tplot, labels)
        plt.tick_params(axis='x', rotation=45, pad=20)
        plt.title('Observation resolution for roi: ' + self.roikeyl[roi] + '(' + titleY + ')')
        plt.xlabel('et')
        plt.ylabel('Resolution [km/px]')
