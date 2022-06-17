from cProfile import label
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, animation
from scipy.fft import fft, ifft


class x_grid:
    def __init__(self, start = 0, end = 500, grid_points = 1000):
        self.grid_points = 1000
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

class potential:
    def __init__(self, x_grid):
        self.x = x_grid.x_values

    def V_values(self):
        
        a = 0.5
        b = 0.01

        v = a*np.exp(-b*(self.x- 150)**2)

        return v

class kinetic_energy:
    def __init__(self, p_grid):
        self.p = p_grid.p_values()

    def KE_values(self):
        ke = (1/2)*self.p**2

        return ke


class td_h:
    def __init__(self, x_grid, potential, kinetic_energy):
        self.x = x_grid
        self.dx = self.x.dx
        self.dt = 0.1
        self.v = potential.V_values()
        self.ke = kinetic_energy.KE_values()

        self.Psi = np.zeros_like(self.x)
    
    def initiallize_Psi(self):
        x = self.x.x_values

        # parameters
        sigma = 5
        x0 = 100
        k = 1

        # constant term 
        c = 1/(sigma*np.sqrt(2*np.pi))

        # mid term
        m = np.exp((-1/2)*((x-x0)/sigma)**2)
        
        # last term
        l = np.exp(1j*k*x)

        self.Psi = c * m * l
        
        return None

    def electric_field(self, N):
        epsilon_0 = 1
        omege = 0.5

        e = epsilon_0 * self.x.x_values * np.cos(omege*N*self.dt)

        return e

    def V_operator(self):
        return np.exp(-0.5 * self.v * self.dt * 1j)
    
    def K_operator(self):
         return np.exp(- self.ke * self.dt * 1j)

    def e_opeartor(self, N):
        epsilon_0 = 3
        omege = 0.5

        e_op = np.exp((-1j*epsilon_0/omege) * 0.5 * self.x.x_values * np.sin(omege * N * self.dt))

        self.Psi *= e_op

        return None

    def split_op(self, N):

        self.Psi *= self.V_operator()

        # FFT to momentum space
        self.Psi = fft(self.Psi)

        # Full step in momentum space
        self.Psi *= self.K_operator()

        # iFFT back
        self.Psi = ifft(self.Psi)

        # Final half-step in real space
        self.Psi *= self.V_operator()



x = x_grid()
p = p_grid(x)
pot = potential(x)
ke = kinetic_energy(p)

tdh = td_h(x, pot, ke)

tdh.initiallize_Psi()



fig, ax = plt.subplots()
plt.close()


ax.set_xlim(( 0, 500))
ax.set_ylim((-0.01, 0.01))
ax.set_xlabel("x values")

line, = ax.plot([], [], label = "Probability Density")
lineV, = ax.plot([], [], label = "Guassian Potential")
ax.legend()

def init():
    line.set_data([],[])
    lineV.set_data([],[])
    return (line, lineV,)

def animate(i):
    tdh.e_opeartor(i)
    tdh.split_op(i)
    tdh.e_opeartor(i)
    ax.set_title("omega = 0.5 and epsilon = 3, t = {time:.2f}".format(time = i*0.1))
    line.set_data(x.x_values, np.abs(tdh.Psi)**2)
    lineV.set_data(x.x_values, pot.V_values())
    return(line,)


anim = animation.FuncAnimation(fig, animate, init_func=init, frames = np.arange(1, 1500),interval=10, blit=True)

anim.save("TDH_11.gif")
