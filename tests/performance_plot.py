import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick


timings = [
    0.04170293807983398,
    0.15635199546813966,
    0.22000808715820314,
    0.3506165027618408,
    0.5720492362976074,
    1.0751613140106202,
    4.725947093963623,
    7.02956256866455,
    103.4705885887146,
]
ncombs = [
    10000,
    160000,
    810000,
    2560000,
    6250000,
    12960000,
    65610000,
    100000000,
    1600000000,
]

fig, axes = plt.subplots(1, 1, figsize=(5, 5))
axes.scatter(ncombs, timings, color="tab:red")
axes.plot(ncombs, timings, color="tab:blue")
# axes.xaxis.set_major_formatter(mtick.FormatStrFormatter("%.1e"))
axes.set_xlim(min(ncombs) * 0.95, max(ncombs) * 1.05)
axes.set_ylim(-1, max(timings) * 1.05)
axes.set_xscale("log")
axes.set_xlabel("log(Number of grid points)")
axes.set_ylabel("Execution time in seconds")
axes.set_title("Scaling with respect to $N_{max}$")

plt.show()
