import numpy as np

def velocity_verlet(type_of_potential, particle, number_of_steps, total_time):
    """This function takes type of potential, particle, number of steps and 
    total time as arguments and return as tuple of x and p values. """

    t_delta = total_time/number_of_steps

    lis_x = []
    lis_p = []

    x_t = particle.get_initial_position()
    v_t = particle.get_initial_velocity()
    p_t = v_t*particle.get_mass()

    lis_x.append(x_t)
    lis_p.append(p_t)
    

    for t in np.arange(0, total_time, t_delta):
        x_t_plusdelta = x_t + t_delta*v_t + (t_delta**2/2)*type_of_potential.get_acc(x_t)
        v_t_plusdelta = v_t + (t_delta/2)*(type_of_potential.get_acc(x_t) + \
            type_of_potential.get_acc(x_t_plusdelta))
        p_t_plusdelta = v_t_plusdelta*particle.get_mass()

        lis_x.append(x_t_plusdelta)
        lis_p.append(p_t_plusdelta)

        x_t = x_t_plusdelta
        v_t = v_t_plusdelta


    return (lis_x, lis_p)

def velocity_verlet_for_OF(type_of_potential, particle, number_of_steps, total_time):
    """This velocity verlet is for oscillating field."""

    t_delta = total_time/number_of_steps

    lis_x = []
    lis_p = []

    x_t = particle.get_initial_position()
    v_t = particle.get_initial_velocity()
    p_t = v_t*particle.get_mass()

    lis_x.append(x_t)
    lis_p.append(p_t)
    

    for t in np.arange(0, total_time, t_delta):

        x_t_plusdelta = x_t + t_delta*v_t + (t_delta**2/2)*type_of_potential.get_acc(t, x_t)
        v_t_plusdelta = v_t + (t_delta/2)*(type_of_potential.get_acc(t, x_t) + \
            type_of_potential.get_acc(t, x_t_plusdelta))
        p_t_plusdelta = v_t_plusdelta*particle.get_mass()

        lis_x.append(x_t_plusdelta)
        lis_p.append(p_t_plusdelta)

        x_t = x_t_plusdelta
        v_t = v_t_plusdelta


    return (lis_x, lis_p)