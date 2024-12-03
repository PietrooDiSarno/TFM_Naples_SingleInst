import copy
import random

import numpy as npy

from FuturePackage.oplanClass import oplan


class oplanSubproblem(oplan):  # a subproblem where only a subset of the dofs of the original problem are changed

    # p: sample of original class;
    # nd: dimensions of the subproblem if it is an integer OR list with the specific dofs to change if a list
    def __init__(self, p, nd):
        if isinstance(nd, int):
            self.nd = nd  # dimensions of the reduced problem
            dl = list(range(p.getNdof()))
            random.shuffle(dl)
            #print('dl=',dl)
            self.newd = dl[0:nd]  # self.newd positions of the original problem to be changed
        elif isinstance(nd, list): # list of the dofs to change
            self.newd = nd
            self.nd = len(self.newd)
        numdim=self.nd
        oplan.__init__(self, self.newd[0], self.newd[-1]) # our unknown now will be a vector with nd components, that will be added to some of the original unknowns
        self.p=p # original class sample
        self.totalDim = p.getNdof() # dimensions of the original problem

    def getChanged(self): # returns a copy of the original class sample with the changes done
        obsLength, stol, qroi, subproblem = self.p.getVector()  # original p converted to a vector
        delta = npy.zeros(self.totalDim)  # delta
        # print('self.newd=',self.newd)
        # print('delta1=',delta)
        # print('x1=',x)
        for i, d in enumerate(self.newd):
            # print(i,d)
            delta[d] = self.stol[i]
        # print('delta2=',delta)
        x = x + delta  # changed vector
        changed = copy.deepcopy(self.p)  # copy of my original
        changed.replaceWithVector(x)  # changed version
        return changed
    def fitFun(self):
        zz=self.getChanged()
        return zz.fitFun()