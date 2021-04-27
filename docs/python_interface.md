# Python Interface

The algorithm can be directly executed from python. This is useful to incorporate the algorithm into other workflows.

```python
from hetbuilder.algorithm import CoincidenceAlgorithm
from hetbuilder.plotting import InteractivePlot

# we read in the structure files via the ASE
import ase.io
bottom = ase.io.read("lower_layer.xyz")
top = ase.io.read("upper_layer.xyz")

# we set up the algorithm class
alg = CoincidenceAlgorithm(bottom, top)
# we run the algorithm for a choice of parameters
results = alg.run(Nmax = 10, Nmin = 0, angles = [0, 10, 15, 20], tolerance = 0.1, weight = 0.5)
```

If the search was not successful, `None` is returned. Otherwise, we can inspect the results like this:
```python
for j in results:
    print(j)
```

This allows to filter the results, e.g., if we only want a certain angle:
```python
k = [j for j in results if j.angle == 0.0]
print(k)
```

We can also parse these results to the matplotlib plotting interface:
```python
iplot = InteractivePlot(bottom=bottom, top=top, results=results, weight=0.5)
iplot.plot_results()
```