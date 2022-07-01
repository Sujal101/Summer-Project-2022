import time

st = time.time()

import numpy as np
import scipy.linalg as sl
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


class tdse:
    def __init__(self) -> None:
        self.x = np.linspace(-50,50,1000)
        self.dt = 0.1
        self.epsilon = 1
        self.omega = 0.1
        self.pot = np.zeros_like(self.x)
        self.iw = np.zeros_like(self.x)
        self.pot_matrix = np.zeros((len(self.x),len(self.x )))
        self.ke_matrix = np.zeros_like(self.pot_matrix)
        self.H_diag = np.zeros_like(self.pot_matrix)
        self.unit_matrix = np.zeros_like(self.H_diag)
        self.inv_unit_matrix =np.zeros_like(self.H_diag)


    def set_potential(self)-> None:
        # DOUBLE-WELL POTENTIAL

        a = 0.003
        b = 0.03
        self.pot = a*self.x**4 - b*self.x**2 +0.05

        #self.pot = 0.5*0.1*self.x**2


        #self.pot = -2.5*np.exp(-0.05*(self.x**2))


        # v = np.zeros_like(self.x)
        # v[500:520] = 0.5
        #self.pot =v        

    def set_initial_wavefunction(self) -> None:
        sigma = 3
        x0 = -10
        k = 1

        # constant term 
        c = 1/(sigma*np.sqrt(2*np.pi))

        # mid term
        m = np.exp((-1/2)*((self.x-x0)/sigma)**2)

        # last term
        l = np.exp(1j*k*self.x)

        self.iw = c * m * l
    
    def set_pot_matrix(self) -> None:
        self.pot_matrix = np.diag(self.pot)

    def set_ke_matrix(self) -> None:
        dx = self.x[1] - self.x[0]
        N = len(self.x)

        cons = 1/(2*(dx**2))
        self.ke_matrix =  np.diag(np.ones(N)*(cons*((np.pi**2)/3)))

        for i in range(N):
            for j in range(N):
                if i != j:
                    self.ke_matrix[i,j] = cons*((-1)**(i-j))*(2/(i-j)**2)
    
    def H0_diag(self):
        h_matrix = self.pot_matrix + self.ke_matrix
        eigv, uni_matrix = sl.eigh(h_matrix)
        h_diag = eigv

        self.H_diag = h_diag
        self.unit_matrix = uni_matrix
        self.inv_unit_matrix = sl.inv(uni_matrix)

    def method(self, N):
        laser_fac =  np.exp(-1j * (self.epsilon/(2*self.omega)) * self.x * np.sin(self.omega * N * self.dt))

        h_fac = np.diag(np.exp(-1j*N*self.dt * self.H_diag))

        wf = laser_fac * self.iw
        wf = self.inv_unit_matrix @ wf.reshape(len(self.x), 1)
        wf = h_fac @ wf
        wf = self.unit_matrix @ wf
        wf = laser_fac * wf.reshape(len(self.x),)

        return wf.reshape(len(self.x),)



t = tdse()
t.set_potential()
t.set_initial_wavefunction()
t.set_pot_matrix()
t.set_ke_matrix()
t.H0_diag()

def animation():
    fig, ax = plt.subplots()

    #ax.set_xlim(t.x[0],t.x[-1])
    ax.set_xlim(t.x[0],t.x[-1])
    ax.set_ylim(-0.05,0.1)
    ax.set_xlabel("x values")

    line, = ax.plot([], [], label= "Probability Distribution of wave packet")
    lineV, = ax.plot([], [], '--', label= "Potential")
    ax.legend()

    def init():
        line.set_data([],[])
        lineV.set_data([],[])

        return (line, lineV,)

    def animate(i):
        fw = t.method(i)

        ax.set_title("omega = {omega: } and epsilon = {epsilon: }, \
            t = {time: .2f}".format(time = i*t.dt, omega = t.omega, epsilon = t.epsilon))
        
        line.set_data(t.x, (np.conjugate(fw)*fw).real)
        lineV.set_data(t.x, t.pot)

        return (line, lineV,)

    anim = FuncAnimation(fig, animate, init_func=init, frames = np.arange(1, 500),interval=10, blit=True)
    anim.save("tdse_dwp4.gif")

#plt.ylim(-0.05,0.1)
# plt.plot(t.x, np.conjugate(t.method(1))*t.method(1))
# plt.show()

animation()
# plt.plot(t.x, np.conj(t.iw)*t.iw)
# plt.plot(t.x, t.pot)
# plt.show()

et = time.time()

print("total time:", et-st)