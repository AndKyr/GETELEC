import gamowCalculator as gc
import time
import matplotlib.pyplot as plt

import getelec_wrap as gt
font = 25
import matplotlib as mb
mb.rcParams["font.size"] = font
mb.rcParams["axes.labelsize"] = font
mb.rcParams["xtick.labelsize"] = font
mb.rcParams["ytick.labelsize"] = font
mb.rcParams["legend.fontsize"] = font
mb.rcParams["lines.linewidth"] = 1.5
mb.rcParams["text.usetex"] = True
figureSize = [16,10]
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

calculator = gc.GamowCalculator(solverType="IVP", XCdataFile="", minimumPotential=15.)

start_time = time.time()
calculator.calculateGamow(4.5)
# calculator.calculateGamowCurve(64, minBarrierDepth=-1., maxGamow=50)
calculator.solver.plotWaveFunction()

emitter = gt.GetelecInterface()