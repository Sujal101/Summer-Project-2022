from particle import particle
from potentials import harmonic_potential, LJ_potential, double_well_potential, morse_potential
from algorithms import velocity_verlet
import matplotlib.pyplot as plt
import random

def phase_space(type_of_potential, particle, number_of_steps, total_time):
    if type_of_potential == 'harmonic':
        val = {"k" : 1.5, "par" : particle}
        pot = harmonic_potential(**val)
        return velocity_verlet(pot, particle, number_of_steps, total_time)

    if type_of_potential == 'morse':
        val = {'d': 100, "a" : 0.04, "par": particle}
        pot = morse_potential(**val)
        return velocity_verlet(pot, particle, number_of_steps, total_time)
    
    if type_of_potential == 'double-well':
        """For this potential to be double-well b and c must be positive numbers."""

        val = {'c': 100, 'b': 1, 'par': particle}
        pot = double_well_potential(**val)
        return velocity_verlet(pot, particle, number_of_steps, total_time)

    if type_of_potential == 'LJ':

        val = {"A": 10, "B":10, "par": particle}
        pot = LJ_potential(**val)
        return velocity_verlet(pot, particle, number_of_steps, total_time)
    
def diff_initial_condi():
    """20 Phase space plots for multiple intial conditions."""

    number_of_plots = 20

    fig, axs = plt.subplots(2,2)
    
    for n in range(number_of_plots):
        initial_position = random.uniform(-10, 10.0)
        initial_positionj = random.uniform(4, 10.0)
        initial_velocity = random.uniform(0,10)

        mass = 5

        par = particle(mass, initial_position, initial_velocity)
        par_j = particle(mass, initial_positionj, initial_velocity)

        (X0, P0) = phase_space("harmonic", par, 10000, 20)
        (X1, P1) = phase_space("morse", par, 10000, 20)
        (X2, P2) = phase_space("double-well", par, 10000, 20)
        (X3, P3) = phase_space("LJ", par_j, 10000, 20)

        axs[0,0].plot(X0,P0)
        axs[0,0].set_xlabel("Position")
        axs[0,0].set_ylabel("Momentum")
        axs[0,0].set_title("Harmonic Potential")
        axs[0,1].plot(X1,P1)
        axs[0,1].set_xlabel("Position")
        axs[0,1].set_ylabel("Momentum")
        axs[0,1].set_title("Morse Potential")
        axs[1,0].plot(X2,P2)
        axs[1,0].set_xlabel("Position")
        axs[1,0].set_ylabel("Momentum")
        axs[1,0].set_title("Double-well Potential")
        axs[1,1].plot(X3,P3)
        axs[1,1].set_xlabel("Position")
        axs[1,1].set_ylabel("Momentum")
        axs[1,1].set_title("LJ Potential")

        
    return plt

diff_initial_condi().show()