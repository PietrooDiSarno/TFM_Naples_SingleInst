import glob
import os
import matplotlib.pyplot as plt
import matplotlib
import numpy as np


from conversion_functions import mat2py_et2utc


class dataHandling:
    def __init__(self, dataPath_ = None):

        if dataPath_:
            self.valPath = os.path.join(dataPath_, 'ROI_makespans', 'makespanValues')
            self.plotPath = os.path.join(dataPath_, 'ROI_makespans', 'plots')
        else:  # Use the path specified with the configuration file or the default
            self.valPath = os.path.join(self.getDefaultDataPath(), 'makespanValues')
            self.plotPath = os.path.join(self.getDefaultDataPath(), 'plots')

        print('dataHandling: Using ', self.valPath, ' for makespan values and ', self.plotPath, ' for plots')

        if not os.path.isdir(self.valPath):
            print('Creating folder', self.valPath)
            os.makedirs(self.valPath, exist_ok = True)

        if not os.path.isdir(self.plotPath):
            print('Creating folder', self.plotPath)
            os.makedirs(self.plotPath, exist_ok = True)

    def getDefaultDataPath(self):
        homeFolder = os.path.expanduser('~')
        dataFolder = os.path.join(homeFolder, 'ROI_makespans')
        return dataFolder

    def getName(self, mosaic, ROIname, n, int):
        return '%s_int%d_%s_%d' % (ROIname, int, mosaic, n)

    def savePlots(self, mosaic, ROIname, start_times, makespan, nImg, int):
        fname = os.path.join(self.plotPath, self.getName(mosaic, ROIname, len(start_times), int) + '.tif')
        plt.figure()
        not_visible = []
        visible = []
        makespan_ = []
        nImg_ = []
        for i,val in enumerate(makespan):
            if val == None:
                not_visible.append(start_times[i])
            else:
                visible.append(start_times[i])
                makespan_.append(makespan[i])
                nImg_.append(nImg[i])

        matplotlib.use('TkAgg')
        fig, ax1 = plt.subplots()
        ax1.set_xlabel('Initial observation instant')
        ax1.set_ylabel('Observation make-span [s]', color='tab:blue')
        ax1.plot(visible, makespan_, 'o', color='b', linestyle='none', markersize=0.75)
        ax1.plot(not_visible, np.zeros(len(not_visible)), 'o', color='r', linestyle='none', markersize=0.75)
        ax1.tick_params(axis='y', labelcolor='tab:blue')
        ax2 = ax1.twinx()
        ax2.set_ylabel('Number of images', color='tab:orange')
        ax2.plot(visible, nImg_, 'o', color='tab:orange', linestyle = 'none', markersize = 0.75)
        ax2.tick_params(axis='y', labelcolor='tab:orange')
        fig.suptitle(f"Makespan and number of images for the {ROIname} ROI. Compliant interval number {int}")
        ax1.set_xticks(np.linspace(start_times[0],start_times[-1],num=6),
                   mat2py_et2utc(np.linspace(start_times[0],start_times[-1],num=6),'C',0), rotation=15)
        fig.tight_layout()
        plt.subplots_adjust(left=0.1, right=0.9, bottom=0.20, top=0.92)
        plt.savefig(fname, format="tiff", bbox_inches = "tight")
        print(f"Plot saved as {self.getName(mosaic, ROIname, len(start_times), int) + '.tif'}")

    def saveValues(self, mosaic, ROIname, start_times, makespan, nImg, int):
        fname = os.path.join(self.valPath, self.getName(mosaic, ROIname, len(start_times), int) + '.txt')
        with open(fname, 'w') as file:
            file.write("start_times:\n")
            file.write(" ".join(map(str, start_times)) + "\n")

            file.write("makespan:\n")
            file.write(" ".join(map(str, makespan)) + "\n")

            file.write("nImg:\n")
            file.write(" ".join(map(str, nImg)) + "\n")

            file.write("num_points:\n")
            file.write(str(len(start_times)) + "\n")
        print(f"File saved as {self.getName(mosaic, ROIname, len(start_times), int) + '.txt'}")

    def readValues(self, filename, ROIname):
        with (open(filename, 'r') as file):
            lines = file.readlines()
            data = {}
            current_key = None
            for line in lines:
                line = line.strip()
                if line.endswith(':'):
                    current_key = line[:-1]
                    data[current_key] = []
                elif current_key:
                    data[current_key].append(line)

            start_times = np.array([float(x) for x in data["start_times"][0].split()])
            makespan = []
            for x in data["makespan"][0].split():
                if x != 'None':
                    makespan.append(float(x))
                else:
                    makespan.append(None)
            makespan = np.array(makespan)

            nImg = []
            for x in data["nImg"][0].split():
                if x != 'None':
                    nImg.append(float(x))
                else:
                    nImg.append(None)
            nImg = np.array(makespan)

            num_points = int(data["num_points"][0])

            print('\n', ROIname, 'ROI','with',num_points,'points','\n', 'start_times:', start_times, '\n','makespan:',
                  makespan, '\n', 'nImg:', nImg, '\n')
            return start_times, makespan, nImg, num_points

    def getValues(self, mosaic, ROIname, int, n_points = None):
        if n_points:
            fname = os.path.join(self.valPath, self.getName(mosaic, ROIname, n_points, int) + '.txt')
            start_times, makespan, nImg, num_points = self.readValues(fname,ROIname)
            return start_times, makespan, nImg, num_points
        else:
            pattern = os.path.join(self.valPath,'%s_int%d_%s' % (ROIname, int, mosaic) + '*.txt')
            file_list = glob.glob(pattern)
            start_times = []
            makespan = []
            nImg = []
            num_points = []
            for fname in file_list:
                aux1, aux2, aux3, aux4 = self.readValues(fname, ROIname)
                start_times.append(aux1)
                makespan.append(aux2)
                nImg.append(aux3)
                num_points.append(aux4)
            return start_times, makespan, nImg, num_points

