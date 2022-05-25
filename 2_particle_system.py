from particle import particle
from algorithms import velocity_verlet
import numpy as np
import matplotlib.pyplot as plt
import random
        
# Two particle in a Harmonic potential and interaction between them is given by Lennard-Jones potential.
# Two particle system reduced to one partice system.
# Here particle mass is reduced mass.

class HP_LJ:

    def __init__(self, **kwargs):
        self.omega = kwargs['omega']
        self.A = kwargs['A']
        self.B = kwargs['B']
        self.par = kwargs['par']
    
    def get_acc(self, x_t):
        a = -(self.omega**2)*x_t - (1/self.par.get_mass())*((self.B/(x_t**7)) \
            - (self.A/(x_t**13)))

        return a


def phase_space(type_of_potential, p1, number_of_steps, total_time):
    if type_of_potential == 'HP_LJ':
        val = {"omega": 1, "A": 10, "B":10, "par": p1}
        pot = HP_LJ(**val)
        return velocity_verlet(pot, p1, number_of_steps, total_time)


def diff_initial_condi(type_of_potential):
    """20 Phase space plots for multiple intial conditions."""

    number_of_plots = 20

    fig, ax = plt.subplots()
    ax.set_xlabel("Position")
    ax.set_ylabel("Momentum")
    ax.set_title(f"Phase space of {type_of_potential} potential with different initial conditions")

    for n in range(number_of_plots):
        initial_position = random.uniform(4, 10.0)
        initial_velocity = random.uniform(0,10)

        p1 = particle(mass=5,initial_position=initial_position, initial_velocity=initial_velocity)

        (x, p) = phase_space(f"{type_of_potential}", p1, 10000, 5)
        plt.plot(x,p)

    return plt


diff_initial_condi("HP_LJ").show()
