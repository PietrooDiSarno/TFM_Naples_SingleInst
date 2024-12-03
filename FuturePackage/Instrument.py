class Instrument:
    def __init__(self, ifov, npix, imageRate, fs):
        self.ifov = ifov
        self.npix = npix
        self.imageRate = imageRate
        self.safetyFactor = fs