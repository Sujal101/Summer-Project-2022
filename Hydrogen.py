import matplotlib.pyplot as plt 
import numpy as np
import random
from matplotlib.animation import FuncAnimation

ALPHA = 5

# I am going to define a function which is eqaul to dr/dt

def dr_dt(r,p):

    dr_dt = p - (r**2)*(p**3)*np.exp(ALPHA*(1-(r*p)**4))

    return dr_dt

# Now a fuction for dp/dt

def dp_dt(r,p):

    dp_dt = -1/r**2 + (1/(2*ALPHA*r**3) + r*p**4)*np.exp(ALPHA*(1-(r*p)**4))

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
     
    #fig, (ax1,ax2,ax3) = plt.subplots(1,3)
    fig, ax3 = plt.subplots()

    # ax1.set_xlabel("time (t)")
    # ax1.set_ylabel("distance between nucleus and electron (r)")
    # ax1.set_title("r vs t")

    # ax2.set_xlabel("time (t)")
    # ax2.set_ylabel("momentum of electron (p)")
    # ax2.set_title("p vs t")

    ax3.set_xlabel("distance between nucleus and electron (r)")
    ax3.set_ylabel("momentum of electron (p)")
    ax3.set_title("p vs r")

    n = 0
    while n < number_of_plots:

        r_initial = random.uniform(1,7)
        p_initial = random.uniform(0,2)
        
        r,p,t = runge_kutta(r_initial, p_initial)

        # ax1.plot(t,r)
        # ax2.plot(t,p)
        ax3.plot(r,p, label="r_i: {r_initial:.2f}, p_i: {p_initial:.2f}".format(r_initial = r_initial, p_initial = p_initial))


        n = n+1

    ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.show()

#with_diff_ini_cond(40)



def animation(number_of_plots, variable):
    n = number_of_plots
    if variable == 'p':
        fig, ax = plt.subplots(constrained_layout=True)
        fig.suptitle('Initial distance between nucleus and electron (r) < 1', fontsize=16)
        ax.set_xlabel("r")
        ax.set_ylabel("p")

        line_lis = []
        r_lis =[]
        for i in range(n):
            r = random.uniform(0,1)
            r_lis.append(r)
            ax.set_xlim(0,50)
            ax.set_ylim(0,5)
            line, = ax.plot([], [], label= "r_initial={x_i:.2f}".format(x_i= r))
            ax.legend()
            line_lis.append(line)

        
        def with_p(p):
            x_lis = []
            y_lis = []
            ax.set_title("p_intial={p:.2f}".format(p=p))
            for i in range(n):
                x,y = runge_kutta(r_lis[i],p)[0:2]
                x_lis.append(x)
                y_lis.append(y)

                line_lis[i].set_data(x_lis[i], y_lis[i])

            return line_lis,
        
        animation = FuncAnimation(fig, with_p, np.arange(0, 2, 0.01))
    
    elif variable == 'r':
        fig, ax = plt.subplots(constrained_layout=True)
        fig.suptitle('Initial momentum of electron (p) < 1', fontsize=16)
        ax.set_xlabel("r")
        ax.set_ylabel("p")

        line_lis = []
        p_lis =[]
        for i in range(n):
            p = random.uniform(0,1)
            p_lis.append(p)
            ax.set_xlim(0,50)
            ax.set_ylim(0,5)
            line, = ax.plot([], [], label= "p_initial={p_i:.2f}".format(p_i= p))
        
            ax.legend()
            line_lis.append(line)

        
        def with_x(r):
            x_lis = []
            y_lis = []
            ax.set_title("r_intial={r:.2f}".format(r=r))
            for i in range(n):
                x,y = runge_kutta(r,p_lis[i])[0:2]
                x_lis.append(x)
                y_lis.append(y)

                line_lis[i].set_data(x_lis[i], y_lis[i])

            return line_lis,
        
        animation = FuncAnimation(fig, with_x, np.arange(1, 2, 0.01))


    plt.show()
    #animation.save("r_less_than_1.gif")

animation(10, 'r')





