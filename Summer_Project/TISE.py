import numpy as np
import scipy.linalg as sl
import matplotlib.pyplot as plt

class grid:

    def __init__(self, start, end, n_points) -> None:
        self.x_initial = start
        self.x_final = end
        self.N = n_points
        self.dx = (end - start)/n_points

    def x_values(self):
        x = np.linspace(self.x_initial, self.x_final, self.N + 1)
        return x[1:-1]

class for_gaussian:

    def __init__(self, grid) -> None:
        self.x = grid.x_values()
        self.dx = grid.dx
        self.N = grid.N
    
    def V_DVR(self):
        v = -0.64*(np.exp(-0.1424*(self.x**2)))

        return np.diag(v)
    
    def T_DVR(self):
        cons = 1/(2*(self.dx**2))
        T_diag =  np.diag(np.ones(self.N-1)*(cons*((np.pi**2)/3)))
        
        for i in range(self.N-1):
            for j in range(self.N-1):
                if i != j:
                    T_diag[i,j] = cons*((-1)**(i-j))*(2/(i-j)**2)

        return T_diag

    def H_DVR(self):

        return self.T_DVR() + self.V_DVR()

class for_mrorse:

    def __init__(self, grid) -> None:
        self.x = grid.x_values()
        self.dx = grid.dx
        self.N = grid.N
    
    def V_DVR(self):

        d = 0.176
        a = 1.04
        xe = 1.4

        v = d*((1-np.exp(-a*(self.x-xe)))**2)

        return np.diag(v)
    
    def T_DVR(self):
        cons = 1/(2*(self.dx**2))
        T =  np.zeros((self.N-1, self.N-1))
        
        for i in range(self.N-1):
            for j in range(self.N-1):
                if i == j:
                    T[i,j] = cons*((np.pi**2)/3 - 1/(2*((i+1)**2)))
                else:
                    T[i,j] = cons*((-1)**(i-j))*((2/(i-j)**2) - (2/(i+j)**2))
        
        return T

    def H_DVR(self):

        return self.T_DVR() + self.V_DVR()


class eigenval_and_eigenvec:

    def __init__(self, h_matrix) -> None:
        self.H = h_matrix
        self.eigenvalues, self.eigenvector_matrix = sl.eigh(h_matrix)
    
    def get_eigenvalues(self):
        return self.eigenvalues
    
    def get_eigenvactor_matrix(self):
        return self.eigenvector_matrix
    
    def get_wavefunction(self, state):
        return self.get_eigenvactor_matrix()[:,state]


def eigenvalues(type_of_potential):

    if type_of_potential == 'morse':
        g = grid(start = 1, end = 100, n_points = 1000)
        x = g.x_values()
        h = for_mrorse(g).H_DVR()
        eigval = eigenval_and_eigenvec(h).get_eigenvalues()

        return eigval
    
    if type_of_potential == 'gaussian':
        g = grid(start = -50, end = 50, n_points = 1000)
        x = g.x_values()
        h = for_gaussian(g).H_DVR()
        eigval = eigenval_and_eigenvec(h).get_eigenvalues()

        return eigval

def plot(state, type_of_potential):

    if type_of_potential == 'morse':
        g = grid(start = 1, end = 100, n_points = 1000)
        x = g.x_values()
        h = for_mrorse(g).H_DVR()

        wave_fun = eigenval_and_eigenvec(h).get_wavefunction(state)

        plt.plot(x, wave_fun)

        return None
    
    elif type_of_potential == 'gaussian':
        g = grid(start = -50, end = 50, n_points = 1000)
        x = g.x_values()
        h = for_gaussian(g).H_DVR()

        wave_fun = eigenval_and_eigenvec(h).get_wavefunction(state)

        plt.plot(x, wave_fun)

        return None
  
plt.plot(0, 'gaussian')
plt.show()