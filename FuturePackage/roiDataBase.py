# This file reads and stores the information regarding the ROIs for JUICE on Callisto and Ganymede
from FuturePackage.roiClass import oPlanRoi
import numpy as np

class ROIDataBase:
    def __init__(self, txt_files=None, bodies=None, customROIs = None):

        #if txt_files is not None and bodies is not None:
        self._ROIs = self._createFromTextFiles(txt_files, bodies)
        #else:
            #self._ROIs = self.createCustomROI(customROIs)

        self._indices = self._getIndices()


    #def createCustomROI(self, customROIs):
    #    myList = []
    #    if not isinstance(customROIs, list):
    #        customROIs = [customROIs]
    #    for ROI in customROIs:
    #        roi['body'] =


    def _createFromTextFiles(self, txt, bodies):
        data = []
        if not isinstance(txt, list):
            txt_files = [txt]
        if not isinstance(bodies, list):
            bodies = [bodies]
        for file, body in zip(txt_files, bodies):
            data += self._parseROIRawData(file, body)
        return data

    def _parseROIRawData(self, file, body):
        mylist = self._readData(file, body)
        cleanlist = self._cleanData(mylist)
        return cleanlist

    def _getIndices(self):
        indices = dict()
        for i, ROI in enumerate(self._ROIs):
            key = ROI['#roi_key']
            indices[key] = i
        return indices

    def getROIs(self, desiredROIs=None):
        myset = set()
        if desiredROIs is None or len(desiredROIs) == 0:
            desiredROIs = self._ROIs
        else:
            desiredROIs = [self._ROIs[self._indices[desiredROIs]]]
        rois = []
        for ROI in desiredROIs:
            if ROI['#roi_key'] in myset:
                print('CAUTION: ROI: ' + ROI + ' has been retrieved twice for the scheduling. Make sure this is '
                                                   'intentional.')
            else:
                myset.add(ROI['#roi_key'])
            i = self._indices[ROI['#roi_key']]
            rois.append(oPlanRoi(self._ROIs[i]['body'],self._ROIs[i]['#roi_key'], self._ROIs[i]['vertices']))
        #if len(rois) == 1:
        #    return rois[0]
        #else:
        return rois
    def getnames(self):
        names = []
        for roi in self._ROIs:
            names.append(roi['#roi_key'])
        return names
    @staticmethod
    def _readData(file, body):
        myset = set()  # To avoid introducing repeated ROIs if any exist on the txt files
        mylist = []
        with open(file, 'r') as myfile:
            header = myfile.readline().strip().split(',')
            for i, name in enumerate(header):
                if name == 'roi_latitudes':
                    header[i] = 'lat'
                elif name == 'roi_longitudes_east':
                    header[i] = 'lon'
            for line in myfile:
                # print(line)
                mydict = dict(zip(header, line.strip().split(',')))
                if mydict['#roi_key'] not in myset:
                    myset.add(mydict['#roi_key'])
                    mydict['body'] = body
                    mylist.append(mydict)
                else:
                    print("ROI: " + mydict['#roi_key'] + ' is repeated on file: ' + file + '. It has been omitted to '
                                                                                           'avoid DataBase repetition')
        return mylist

    @staticmethod
    def _cleanData(mylist):
        headers = ['lat', 'lon']
        for i, mydict in enumerate(mylist):
            for name in headers:
                aux = mydict[name]
                cleanData = aux.strip('[]').split()
                coords = [float(x) for x in cleanData]
                if name == 'lon':
                    for j, coord in enumerate(coords):
                        if coord>180.:
                            coords[j] = coord-360.
                mydict[name] = coords
            mydict['vertices'] = np.array([list(coord) for coord in zip(mydict['lon'], mydict['lat'])])
            mylist[i] = mydict
        return mylist


