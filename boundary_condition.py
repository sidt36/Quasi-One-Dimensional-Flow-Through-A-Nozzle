import numpy as np

def endpoints(density,temperature,velocity):
    #To set the boundary conditions  
    #Density and Temperature are fixed.
    density[:,0] = 1
    temperature[:,0] = 1
    #We can't really call the entry point a reservoir as no mass will flow if the entry velocity is zero, we have to float velocity 
    velocity[:,0] = 2*velocity[:,1] - velocity[:,2]
    
    #All three flow variables have to float in the exit
    density[:,-1] = 2*density[:,-2] - density[:,-3]
    temperature[:,-1] = 2*temperature[:,-2] - temperature[:,-3]
    velocity[:,-1] = 2*velocity[:,-2] - velocity[:,-3]
    
    return density,temperature,velocity