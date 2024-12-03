import subprocess
import time

class spawnProcess:

    def __init__(self, nroi, runName, nroiprocess = 1):
        self.nroi = nroi
        self.nroiprocess = nroiprocess
        self.runName = runName
        self.nprocess = nroi // nroiprocess
        self.nroirem = nroi % nroiprocess

    def whileRunning(self):
        time.sleep(4)
        print('-----whileRunning ')

    def spawn(self):
        subp = []
        finished = []
        if self.nroirem != 0:
            processnum = self.nprocess + 1
        else:
            processnum = self.nprocess
        for process in range(processnum):
            if process <= self.nprocess - 1:
                numroiprocess = self.nroiprocess
            else:
                numroiprocess = self.nroirem

            command = ["python", self.runName + '.py', str(process + 1), str(numroiprocess),str(self.nroirem)]
            print(command)
            subp.append(subprocess.Popen(command))
            #    subp.append(subprocess.Popen(["sleep", str(agent)]))
            finished.append(0)

        while True:
            nf = 0
            for a in range(0, processnum):
                ret = subp[a].poll()
                if ret is not None:
                    finished[a] = 1
                    print('**** process', a + 1, 'finished ', ret)
                    nf = nf + 1
            if nf == processnum:
                print('**** all finished')
                break

            self.whileRunning()


