class particle:
    def __init__(self, mass, initial_position, initial_velocity):
        """It takes three arguments mass of the particle, its initial position 
        and its initial velocity."""

        self.mass = mass
        self.initial_position = initial_position
        self.initial_velocity = initial_velocity
    
    def get_mass(self):
        return self.mass
    
    def get_initial_position(self):
        return self.initial_position
    
    def get_initial_velocity(self):
        return self.initial_velocity