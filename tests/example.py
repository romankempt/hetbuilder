from hetbuilder.algorithm import CoincidenceAlgorithm
from hetbuilder.plotting import InteractivePlot

# we read in the structure files via the ASE
import ase.io
bottom = ase.io.read("graphene.xyz")
top = ase.io.read("MoS2_2H_1l.xyz")

# we set up the algorithm class
alg = CoincidenceAlgorithm(bottom, top)
# we run the algorithm for a choice of parameters
results = alg.run(Nmax = 12, Nmin = 0, angles = [0, 10, 15, 20, 25, 30], tolerance = 0.2, weight= 0.5)

iplot = InteractivePlot(bottom=bottom, top=top, results=results, weight=0.5)
iplot.plot_results()
