from multiprocessing.spawn import import_main_path
from particle import particle
from potentials import electrostatic_potential
from algorithms import velocity_verlet
import matplotlib.pyplot as plt
import numpy as np

par = particle(1, 1, 1)
val = {"Z":1, "particle": par}
pot = electrostatic_potential(**val)

(x,p) = velocity_verlet(pot, par, 10000, 5.75)

plt.plot(x,p)
plt.show()
