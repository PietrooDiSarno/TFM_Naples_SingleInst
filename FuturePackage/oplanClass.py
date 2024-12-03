import random
from pySPICElib.SPICEtools import *
from PSOA.pointres import pointres

from FuturePackage import DataManager


class oplan():
    def __init__(self, substart = None, subend = None):
        self.subproblem = [substart, subend]
        roiL = DataManager.getInstance().getROIList(self.subproblem[0], self.subproblem[1])
        self.stol = [0.] * len(roiL)  # start time observation list
        self.qroi = [0.] * len(roiL)  # cache for the fitness of each ROI observation
        self.obsLength = [0.] * len(roiL)  # cache for the duration of the observation for each ROI

    def removeEmptyTW(self, unconstrainedTW):
        roiL = DataManager.getInstance().getROIList(self.subproblem[0], self.subproblem[1])
        tws = []
        if spice.wncard(unconstrainedTW) == 1:
            print("Unconstrained Search Space is a unique TW")
            return unconstrainedTW
        for tw in range(spice.wncard(unconstrainedTW)):
            twbeg, twend = self.getTWBeginEnd(tw)
            isempty = True
            for roi in roiL:
                if not isempty: break
                for interval in range(spice.wncard(roi)):
                    twStart, twEnd = spice.wnfetd(roi, interval)
                    if twStart >= twbeg and twEnd <= twend:
                        isempty = False
                        tws.append(tw)
                        break
        r = stypes.SPICEDOUBLE_CELL(len(tws) * 2)
        for i in range(len(tws)):
            t0, t1 = spice.wnfetd(unconstrainedTW, i)
            spice.wninsd(t0, t1, r)
        return r


    def getTWBeginEnd(self, tw, interval):  # ARREGLAR returns the et begin & end of a TW interval
        TWBegin, TWEnd = spice.wnfetd(tw, interval)
        return TWBegin, TWEnd

    def print_auxdata(self):
        roiL = DataManager.getInstance().getROIList(self.subproblem[0], self.subproblem[1])
        instrument = DataManager.getInstance().getInstrumentData()
        observer = DataManager.getInstance().getObserver()

        print('observer:', observer)
        print('target:', roiL[0].body)
        print('Regions to be studied:', )
        for roi in roiL:
            print('     -', roi.name)
        print('ROIS in this plan ', len(roiL))
        for i, roi in enumerate(roiL):
            print_tw(roi.ROI_TW, '(s) for ' + roi.name)
        print('ifov=', instrument.ifov)

    def getAllTw(self, i):  #   returns a single interval tw plus start, end covering all tw;
        roiL = DataManager.getInstance().getROIList(self.subproblem[0], self.subproblem[1])
        sa = float('inf')  # first start
        ea = -float('inf')  # last end
        TWbeg, TWend = self.getTWBeginEnd(i)
        for tw in roiL:  # NOT LIKE THIS
            for interval in range(spice.wncard(tw)):
                s, e = spice.wnfetd(tw, interval)
                if s >= TWbeg and e <= TWend:  # If the interval of the general timewindow (the one that contains info for all fbs and rois) is within the flyby we are studying
                    if s < sa: sa = s
                    if e > ea: ea = e
                if sa != float('inf') and ea != -float('inf'):
                    r = newTimeWindow(sa, ea)
        return r, sa, ea

    def getObsLength(self, roi, et):
        interval, _, _ = self.findIntervalInTw(et, roi.ROI_TW)
        _, timeobs, _ = roi.interpolateObservationData(et, interval)
        return timeobs

    # evals a metric of Res of all ROIs (the LOWER, the better)
    # in this case, averaged centroid res. in km/pix

    def evalResPlan(self):
        for i in range(len(self.stol)):
            ts = self.stol[i]
            te = ts + self.obsLength[i]
            et = np.linspace(ts, te, 4)
            qv = []
            for t in et:
                qv.append(self.evalResRoi(i, t))
            self.qroi[i] = sum(qv) / len(et)

        return self.qroi

    def evalResRoi(self, i, et):  # returns instantaneous resolution (fitness) of roi (integer)
        roiL = DataManager.getInstance().getROIList(self.subproblem[0], self.subproblem[1])
        observer = DataManager.getInstance().getObserver()
        instrument = DataManager.getInstance().getInstrumentData()
        print(i)
        _, _, res = roiL[i].interpolateObservationData(et)
        return res  # pointres(instrument.ifov, roiL[i].centroid, et, roiL[i].body, observer)

    # returns the total overlap time, defined as the sum of the overlaps between consecutive observations
    # overlaps of more that two observations are not considered
    def getTotalOverlapTime(self):
        si = sorted(range(len(self.stol)), key=lambda k: self.stol[k])  # ROI sorted by obs time
        sorted_stol = [self.stol[i] for i in si]
        isfirst = True
        toverlap = 0
        for i in range(len(sorted_stol)):
            if isfirst:
                isfirst = False
                continue
            startt = sorted_stol[i]
            endprevious = sorted_stol[i - 1] + self.obsLength[si[i - 1]]
            overlap = 0
            if endprevious > startt: overlap = endprevious - startt
            toverlap = toverlap + overlap
        return toverlap

    def ranFun(self):
        roiL = DataManager.getInstance().getROIList(self.subproblem[0], self.subproblem[1])
        for i, roi in enumerate(roiL):
            _, rr, obslen = self.uniformRandomInTw(roi)
            self.stol[i] = rr
            self.obsLength[i] = obslen

    def mutFun(self, f=0, g=0):
        roiL = DataManager.getInstance().getROIList(self.subproblem[0], self.subproblem[1])
        for i, roi in enumerate(roiL):
            currentBegin = self.stol[i]
            newBeginEt, obslen = self.randomSmallChangeIntw(currentBegin, roi, f)
            self.stol[i] = newBeginEt
            self.obsLength[i] = obslen

    def fitFun(self):
        tov = self.getTotalOverlapTime()
        if tov > 0:
            return tov * 1e9
        return np.mean(self.evalResPlan())


    def uniformRandomInTw(self, roi):
        nint = spice.wncard(roi.ROI_TW)
        plen = [0] * nint
        outOfTW = True

        for i in range(nint):
            intbeg, intend = spice.wnfetd(roi.ROI_TW, i)
            plen[i] = intend - intbeg
        total = sum(plen)
        probabilities = [p / total for p in plen]
        val = list(range(nint))
        # Select one float randomly with probability proportional to its value
        psel = random.choices(val, weights=probabilities)[0]
        i0, i1 = spice.wnfetd(roi.ROI_TW, psel)
        while outOfTW:
            rr = random.uniform(i0, i1)
            obslen = self.getObsLength(roi, rr)
            if rr + obslen <= i1:
                outOfTW = False
        return psel, rr, obslen

    #     # given time t and a tw with multiple intervals, returns the number of interval contaning t, the start and the end
    #     # returns -1,0.0,0.0 if t is not contained in any interval
    def findIntervalInTw(self, t, tw):
        nint = spice.wncard(tw)
        for i in range(nint):
            intbeg, intend = spice.wnfetd(tw, i)
            if intbeg <= t <= intend:
                return i, intbeg, intend
        return -1, 0.0, 0.0

    # given time t0, with both t0 and t0+olen belonging to the same interval in time window tw,
    # returns another instant t0new so that t0new and t0new+olen are in the same tw interval as before
    # t0new-t0 distributes N(0,sigma), when the conditions are satisfied

    def randomSmallChangeIntw(self, t0, roi, f):
        i, intervals, intervalend = self.findIntervalInTw(t0, roi.ROI_TW)
        if np.abs(t0 - intervals) >= np.abs(t0 - intervalend):
            sigma0 = np.abs(intervalend - t0)
        else:
            sigma0 = np.abs(intervals - t0)
        sigma = sigma0
        ns = 0
        while True:
            newBegin = t0 + np.random.normal(0, sigma)
            obslen = self.getObsLength(roi, newBegin)
            newEnd = newBegin + obslen
            # print('NewEnd ', newEnd)
            # print('Intervalend', intervalend)
            if newBegin >= intervals and newEnd <= intervalend:
                # print('Mutation with sigma = ' + str(sigma) + 's')
                break
            ns += 1
            # print('iteration ', ns)
            if ns > 50:
                # print('halving')
                sigma0 = sigma0 / 2
                sigma = sigma0
            if ns > 500:
                print(t0)
                print(intervals)
                print(intervalend)
                raise Exception('uhhh cant find mutation')

        return newBegin, obslen

    def distance(self, other):
        dd = 0
        for i in range(len(self.stol)):
            q = math.fabs(other.stol[i] - self.stol[i])
            if q > dd: dd = q
        return q


    def repFun(self, p1, f1, f2):
        roiL = DataManager.getInstance().getROIList(self.subproblem[0], self.subproblem[1])
        newind = []
        newobs = []
        for i in range(len(p1.stol)):
            op = (p1.stol[i] + self.stol[i])/2
            #print(len(p1.stol))
            #print(len(self.stol))
            a,_,_ = self.findIntervalInTw(op, roiL[i].ROI_TW)
            if a != -1:
                #print('PARENT FOUND')
                newind.append(op)
                obslen = self.getObsLength(roiL[i], op)
                newobs.append(obslen)
            else:
                #print('Op not found, mutating')
                _, rr, obslen = self.uniformRandomInTw(roiL[i])
                newind.append(rr)
                newobs.append(obslen)
        self.stol = newind
        self.obsLength = newobs




    def getNdof(self):
        return len(self.stol)

    def getVector(self):
        return self.obsLength, self.stol, self.qroi, self.subproblem

    def replaceWithVector(self, newvect):
        self.stol = newvect
