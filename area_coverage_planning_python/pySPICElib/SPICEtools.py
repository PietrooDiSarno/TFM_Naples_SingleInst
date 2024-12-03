import numpy as np

import spiceypy as spice
import spiceypy.utils.support_types as stypes
import matplotlib.pyplot as plt
import matplotlib.patches as patches


def print_tw(tw, label=None, detailed=False):
    nint = spice.wncard(tw)
    if label:
        print('time window', label, end=' ')
    print('number of intervals', nint)
    for i in range(nint):
        intbeg, intend = spice.wnfetd(tw, i)
        print('interval', i, 'start ', intbeg,'=',spice.et2datetime(intbeg), 'end ', intend,'=',spice.et2datetime(intend), 'length (h)=',
              (intend - intbeg) / 3600.0)

def plot_tw(ax,tw,lo,hi,color):
    nint = spice.wncard(tw)
    for i in range(nint):
        intbeg, intend = spice.wnfetd(tw, i)
        ax.add_patch(patches.Rectangle( (intbeg,lo),width=intend-intbeg, height=hi-lo,lw=1,color=color,fill=True))

def etToAxisStrings(et, ndates=5, accurate=False):
    # Select vector of et
    if isinstance(ndates, int):  # just the number of marks
        idx = np.round(np.linspace(0, len(et) - 1, ndates)).astype(int)
        etv = [et[i] for i in idx]  # we select a ndates equally spaced values
        convert = True
    else:  # a list, we assume they are strings and generate the et
        convert = False
        etv = [spice.str2et(d) for d in ndates]
        ets = ndates

    # Convert to strings (if needed)
    if convert:
        ets = []
        for x in etv:
            ww = spice.et2utc(x, 'C', 0)
            ww = ww.split(' ')
            s = ww[0] + ' ' + ww[1] + ' ' + ww[2] + ' '+ ww[3]
            ets.append(s)

    return etv, ets


# userfun: user defined fuction, takes one float (et) and returns one float
# cnfine: time window to search OR list with two floats that are the start and end of et for the search
# relate: '<' or '>'
# refval: value searched (eg., < 1e6 )
# step: search interval as in SPICE functions (see doc)
# max number of intervals of the tw returned
# returns a time window

def mySolver(userfun,e0,e1,nite,tolerance): # Auxiliary of myTwFinder (Bolzano slow, by now)
    a=e0
    fa=userfun(e0)
    b=e1
    fb=userfun(e1)
    ite=0
    if fa*fb>0:
        raise Exception('fa*fb>0')

    while True:
        c=(a+b)/2.
        fc=userfun(c)
        if fa*fc>0:
            a=c
        else:
            b=c
        if abs(fc)<tolerance or ite>nite:
            return c,fc
        ite=ite+1

def myTwFinder(userfun, cnfine, relate, refval, step, maxntw=200):
    info=0
    def val(x): # val<0: false; >0: true
        q = userfun(x) - refval
        if relate=='<': q = -q
        return q

    r = stypes.SPICEDOUBLE_CELL(maxntw) # tw to contain the result

    if isinstance(cnfine,spice.utils.support_types.SpiceCell):
        pass
    else:
        t0=cnfine[0]
        t1=cnfine[1]
        cnfine = stypes.SPICEDOUBLE_CELL(2)
        spice.wninsd(t0, t1, cnfine)
    nint = spice.wncard(cnfine)
    if info: print('number of intervals in search time window', nint)
    for i in range(nint): # for each interval in the search time window
        et0, et1 = spice.wnfetd(cnfine, i) # get start & end
        rr=np.arange(et0+step,et1+step,step)

        if val(et0)>0:
            if info: print('>> START1: ',et0, val(et0),spice.et2datetime(et1))
            tstart=et0
            state=True
        else:
            state=False
        for t in rr:
            v=val(t)
            if v>0: # is true
                if state==False: #and was false, cioè al tempo et0 val(et0)<0 per la prima iterazione del ciclo for
                    zt, _ = mySolver(val, t - step, t, 20, 0.01)
                    if info: print('>> START2: ', zt,val(zt),spice.et2datetime(zt))
                    tstart = zt
                state = True
            else: # is false
                if state==True: # and was true
                    zt,_ = mySolver(val, t - step, t, 20, 0.01)
                    if info: print('>> END1: ',t,val(t),spice.et2datetime(zt))
                    tend = zt
                    spice.wninsd(tstart, tend, r)  # add an interval
                state = False
        if state==True:
            if info: print('>> END2: ', et1, val(et1),spice.et2datetime(et1))
            tend = et1
            spice.wninsd(tstart, tend, r)

    return r



