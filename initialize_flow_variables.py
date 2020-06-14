import numpy as np
from boundary_condition import endpoints

def initialize_with_zeros(timesteps,spacesteps):
    density = np.zeros((timesteps,spacesteps))
    temperature = np.zeros((timesteps,spacesteps))
    velocity = np.zeros((timesteps,spacesteps))
    # We create a placeholder for our data
    return density,temperature,velocity

def initialize_with_values(delta_x,timesteps,spacesteps):
    (density , temperature , velocity)= initialize_with_zeros(int(timesteps),int(spacesteps)) 
    #We provide linearly varying values to out flow variables at time = 0 to ensure easy convergence
    for i in range(spacesteps):
        density[0,i] = 1 - 0.3146*delta_x*i
        temperature[0,i] = 1 - 0.2314*delta_x*i
        velocity[0,i] = (0.1 + 1.09*delta_x*i)*temperature[0,i]**(0.5)
        
    density,temperature,velocity = endpoints(density,temperature,velocity)
    return density,temperature,velocity                                         
