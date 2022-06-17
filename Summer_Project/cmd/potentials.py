from keyword import kwlist
import numpy as np


class harmonic_potential:
    """form of harmonic potential is v(x) = (1/2)kx^2."""

    def __init__(self, **kwargs):
        self.k = kwargs['k']
        self.par = kwargs['par']
    
    def get_force(self, x_t):
        f = -self.k*x_t
        return f

    def get_acc(self, x_t):
        return self.get_force(x_t)/self.par.get_mass()

class LJ_potential:

    def __init__(self, **kwargs):
        self.A = kwargs["A"]
        self.B = kwargs["B"]
        self.par = kwargs["par"]

    def get_force(self, x_t):
        f = (-(self.B/(x_t**7)) + (self.A/(x_t**13)))
        return f

    def get_acc(self, x_t):
        return self.get_force(x_t)/self.par.get_mass()

class double_well_potential:
    """form of double-well potential is v(x) = bx^4 - cx^2."""

    def __init__(self, **kwargs):
        self.c = kwargs['c']
        self.b = kwargs['b']
        self.par = kwargs['par']

    def get_force(self, x_t):
        f = 2*self.c*x_t - 4*self.b*(x_t)**3
        return f
    
    def get_acc(self, x_t):
        return self.get_force(x_t)/self.par.get_mass()


class morse_potential:
    """form of morse potential is v(x) = d(1-exp(-ax))^2."""

    def __init__(self, **kwargs):
        self.d = kwargs['d']
        self.a = kwargs['a']
        self.par = kwargs['par']
    
    def get_force(self, x_t):
        f = -2*self.a*self.d*(1-np.exp(-self.a*x_t))*(np.exp(-self.a*x_t))
        return f
    
    def get_acc(self, x_t):
        return self.get_force(x_t)/self.par.get_mass()

class electrostatic_potential:
    def __init__(self, **kwargs) -> None:
        self.Z = kwargs["Z"]
        self.par = kwargs["particle"]
        
    def get_force(self, x_t):
        f = -self.Z/x_t**2
        return f
    
    def get_acc(self, x_t):
        return self.get_force(x_t)/self.par.get_mass()

class heisenberg_potential:
    def __init__(self, **kwargs) -> None:
        self.c = kwargs["c"]
        self.b = kwargs["b"]
        self.par = kwargs["particle"]
        

    def get_force(self, x_t, p_t):
        f = (2*self.c/x_t**3 + 4*self.b*(p_t**4) * x_t**3)*np.exp(-self.b*(x_t*p_t)**4)
        return f
    
    def get_acc(self, x_t, p_t):
        return self.get_force(x_t, p_t)/self.par.get_mass()