import matplotlib.pyplot as plt 
import numpy as np
import random

ALPHA = 5

# I am going to define a function which is eqaul to dr/dt

def dr_dt(r,p):

    dr_dt = p - (r**2)*(p**3)*np.exp(ALPHA*(1-(r*p)**4))

    return dr_dt

# Now a fuction for dp/dt

def dp_dt(r,p):

    dp_dt = 1/r**2 + (1/(2*ALPHA*r**3) + r*p**4)*np.exp(ALPHA*(1-(r*p)**4))

    return dp_dt

# Algorithm for intergrating equations of motion

def runge_kutta(initial_position, initial_momentum, starting_time=0, final_time=20, time_step=0.01):
    """This function takes initial position, initial momentum, starting time, final time 
    and time step as arguments and return a list of values of positon, values of momentum 
    and time values."""

    r_lis = []
    p_lis = []
    t_lis = []

    dt = time_step # step size

    r = initial_position
    p = initial_momentum

    t = starting_time

    while t <= final_time:
        k1_r = dr_dt(r,p)
        k1_p = dp_dt(r,p)
        k2_r = dr_dt(r+dt*k1_r/2, p+dt*k1_p/2)
        k2_p = dp_dt(r+dt*k1_r/2, p+dt*k1_p/2)
        k3_r = dr_dt(r+dt*k2_r/2, p+dt*k2_p/2)
        k3_p = dp_dt(r+dt*k2_r/2, p+dt*k2_p/2)
        k4_r = dr_dt(r+dt*k3_r, p+dt*k3_p)
        k4_p = dp_dt(r+dt*k3_r, p+dt*k3_p)

        r = r + dt/6*(k1_r + 2*k2_r + 2*k3_r + k4_r)
        p = p + dt/6*(k1_p + 2*k2_p + 2*k3_p + k4_p)

        r_lis.append(r)
        p_lis.append(p)

        t = t + dt
        t_lis.append(t)

    return [r_lis, p_lis, t_lis]


# plots with different initial conditions

def with_diff_ini_cond(number_of_plots):
     
    fig, (ax1,ax2,ax3) = plt.subplots(1,3)

    ax1.set_xlabel("time (t)")
    ax1.set_ylabel("distance between nucleus and electron (r)")
    ax1.set_title("r vs t")

    ax2.set_xlabel("time (t)")
    ax2.set_ylabel("momentum of electron (p)")
    ax2.set_title("p vs t")

    ax3.set_xlabel("distance between nucleus and electron (r)")
    ax3.set_ylabel("momentum of electron (p)")
    ax3.set_title("p vs r")

    legend_lis_ax1 = []
    legend_lis_ax2 = []
    legend_lis_ax3 = []

    n = 0
    while n < number_of_plots:

        r_initial = random.uniform(-5,0)
        p_initial = random.uniform(-1,1)

        legend_lis_ax1.append("r_initial: {r_initial:.2f}".format(r_initial = r_initial))
        legend_lis_ax2.append("p_initial: {p_initial:.2f}".format(p_initial = p_initial))
        legend_lis_ax3.append("r_initial: {r_initial:.2f}\np_initial: {p_initial:.2f}".format(r_initial = r_initial, p_initial = p_initial))
        
        r,p,t = runge_kutta(r_initial, p_initial)

        ax1.plot(t,r)
        ax2.plot(t,p)
        ax3.plot(r,p)

        n = n+1


    ax1.legend(legend_lis_ax1)
    ax2.legend(legend_lis_ax2)
    ax3.legend(legend_lis_ax3)

    plt.show()

with_diff_ini_cond(5)


# r,p,t = runge_kutta(-5, 1, 0, 10, 0.01)

# plt.plot(t,r)
# plt.plot(t,p)
# plt.plot(r,p)
# plt.legend(['r vs t', 'p vs t', 'p vs r'])
# plt.show()
