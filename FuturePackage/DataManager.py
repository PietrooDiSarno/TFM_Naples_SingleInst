import copy
import numpy as np

# PONER COMENTARIOS
class DataManager:
    __instance = None
    __lock = False

    @staticmethod
    def getInstance():
        if DataManager.__instance is None:
            DataManager()
        return DataManager.__instance

    def __init__(self, roiL, instrumentData, observer):
        if roiL is not None and instrumentData is not None and observer is not None:
            if DataManager.__lock:
                raise Exception("Already initialized")
            self.roiList = roiL
            self.instrument = instrumentData
            self.observer = observer
            self.maxRes = None
            self.minRes = None
            self.getMaxMinRes()
            DataManager.__lock = True
        DataManager.__instance = self


    def getROIList(self,s = None, e = None):
        if s is None:
            return self.roiList
        else:
            return self.roiList[s:e+1]

    def getMaxMinRes(self):
        if self.maxRes == None or self.minRes == None:
            min_res = np.inf
            max_res = - np.inf
            for roi in self.roiList:
                min_res_roi = min(np.concatenate(roi.ROI_ObsRes))
                max_res_roi = max(np.concatenate(roi.ROI_ObsRes))
                if min_res_roi < min_res: min_res = copy.deepcopy(min_res_roi)
                if max_res_roi > max_res: max_res = copy.deepcopy(max_res_roi)
            self.minRes = copy.deepcopy(min_res)
            self.maxRes = copy.deepcopy(max_res)
        else:
            return self.maxRes, self.minRes
    def getSingleROI(self, i):
        return self.roiList[i]

    def getInstrumentData(self):
        return self.instrument

    def getObserver(self):
        return self.observer

    def initObservationDataBase(self, twL, i = None):
        if i is None: #i is the roi number, if None it is assumed that a sorted list containing the constrained TW for each roi is being passed. It must have the same order as the roiList
            for i, tw in enumerate(twL):
                self.roiList[i].initializeObservationDataBase(tw, self.instrument, self.observer, mosaic = True)
        else:
            self.roiList[i].initializeObservationDataBase(twL, self.instrument, self.observer, mosaic = True)