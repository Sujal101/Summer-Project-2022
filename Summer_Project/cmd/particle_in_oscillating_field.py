from particle import particle
from algorithms import velocity_verlet_for_OF
import numpy as np
import matplotlib.pyplot as plt
import random

class oscillating_field:
    """form of oscillating field is \
        v(x,t) = -1/|x| + epsilon_0 * x * cos(omega * t)"""

    def __init__(self, **kwargs) -> None:
        self.epsilon_0 = kwargs["epsilon_0"]
        self.omega = kwargs["omega"]
        self.par = kwargs["particle"]
    
    def get_force(self, t, x_t):
        if x_t > 0:
            f = 1/(x_t)**2 + self.epsilon_0 * np.cos(self.omega * t)
        elif x_t < 0:
            f = - 1/(x_t)**2 + self.epsilon_0 * np.cos(self.omega * t)
        else:
            print("x = 0 can not allowed.")

        return f
    
    def get_acc(self, t, x_t):
        return self.get_force(t, x_t)/self.par.mass


def diff_initial_cond():

    number_of_plots = 10

    fig, axs = plt.subplots(3,3)
    epsilon_0 = 5
    omega = 3

    for j in range(3):
        for k in range(3):
            axs[j,k].set_xlabel("Position")
            axs[j,k].set_ylabel("Momentum")

            epsilon_0 = round(random.uniform(0,10), 2)
            omega = round(random.uniform(0,3), 2)

            axs[j,k].set_title(f"Phase space with different initial conditions\n \
                epsilon_0: {epsilon_0} and omega: {omega}")
            

            for i in range(number_of_plots):
                initial_position = round(random.uniform(-10,10), 2)
                initial_velocity = round(random.uniform(-5, 5), 2)

                p1 = particle(mass=5, initial_position=initial_position, initial_velocity=initial_velocity)
                val = {"epsilon_0": epsilon_0, "omega": omega, "particle": p1}
                pot = oscillating_field(**val)
                (x, y) = velocity_verlet_for_OF(pot, p1, 1000, 5)
                axs[j,k].plot(x,y)

    plt.show()

diff_initial_cond()