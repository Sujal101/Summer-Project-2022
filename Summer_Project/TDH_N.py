import time

st = time.time()


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.fft import fft, ifft



class x_grid:

    def __init__(self, start = 0, end = 500, grid_points = 1000) -> None:
        self.x_values = np.linspace(start, end, grid_points)
        self.dx = self.x_values[1] - self.x_values[0]

class p_grid:
    def __init__(self, x_grid):
        self.x = x_grid.x_values
        self.dx = x_grid.dx

    def p_values(self):
        res = len(self.x)
        dp = 2 * np.pi / (res * self.dx)
        p = np.concatenate((np.arange(0, res / 2),
                                    np.arange(-res / 2, 0))) * dp
        
        return p
    

# t = t0 + N*dt

def initial_wavefunction(x_values):

        # parameters
        sigma = 5
        x0 = 100
        k = 1

        # constant term 
        c = 1/(sigma*np.sqrt(2*np.pi))

        # mid term
        m = np.exp((-1/2)*((x_values-x0)/sigma)**2)
        
        # last term
        l = np.exp(1j*k*x_values)

        Psi = c * m * l
        
        return Psi


def final_wavefunction(wavefunction, x_values, p_values, N, epsilon, omega, dt):

    def TDP_op(x_values, N, epsilon, omega, wavefunction):
        val = np.exp((-1j*epsilon/(2*omega)) * x_values * np.sin(omega * N * dt))

        return val*wavefunction


    def TIP_op(x_values, p_values, N , wavefunction):

        def potential_op(wavefunction):

            def potential():
                v = np.zeros_like(x_values)

                v[300:330] = 0.7

                return v

            val = np.exp((-1j/2)*dt*potential())

            return val * wavefunction

        def kinetic_op(wavefunction):

            def kinetic():

                val =  0.5 * p_values**2

                return val
            
            val = np.exp(-1j * dt * kinetic())

            return val * wavefunction

        for i in range(N):
            wavefunction1 = potential_op(wavefunction)

            # position basis to momentum basis
            wavefunction1 = fft(wavefunction1) 

            wavefunction2 = kinetic_op(wavefunction1)

            # momentum basis to position basis
            wavefunction2 = ifft(wavefunction2) 

            wavefunction3 = potential_op(wavefunction2)

            wavefunction = wavefunction3
        
        return wavefunction

    wavefunction1 = TDP_op(x_values, N, epsilon, omega, wavefunction)
    wavefunction2 = TIP_op(x_values, p_values, N, wavefunction1)
    wavefunction3 = TDP_op(x_values, N, epsilon, omega, wavefunction2)

    return wavefunction3


# Define x and p values

x = x_grid()
x_values = x.x_values
p = p_grid(x)
p_values = p.p_values()
iw = initial_wavefunction(x_values)
def potential_plot(x_values):
    v = np.zeros_like(x_values)

    v[300:330] = 0.7

    return v
pot = potential_plot(x_values)


# NOW LET'S START ANIMATION


fig, ax = plt.subplots()

ax.set_xlim(( 0, 500))
ax.set_ylim((-0.01, 0.01))
ax.set_xlabel("x values")

line, = ax.plot([], [], label= "Probability Distribution")
lineV, = ax.plot([], [], label= "Potential")
ax.legend()

def init():
    line.set_data([],[])
    lineV.set_data([],[])

    return (line, lineV,)

def animate(i):
    epsilon = 3
    omega = 1
    dt = 0.1
    fw = final_wavefunction(iw, x_values, p_values, N = i, \
        epsilon=epsilon, omega=omega, dt = dt)

    ax.set_title("omega = {omega: } and epsilon = {epsilon: }, \
        t = {time:.2f}".format(time = i*dt, omega = omega, epsilon = epsilon))
    
    line.set_data(x_values, np.abs(fw)**2)
    lineV.set_data(x_values, pot)

    return (line, lineV,)

anim = FuncAnimation(fig, animate, init_func=init, frames = np.arange(1, 1000),interval=10, blit=True)
anim.save("TDH_N7.gif")

ed = time.time()

print("total time taken by the code:", ed-st)