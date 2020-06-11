import numpy as np
def speed_non_dimensional(temp_at_x):
    #SI Units
    speed =  (temp_at_x)**0.5
    return speed


def time_step_taken(temperature_at_t,courant_number,deltax,velocity_at_t):
    small = float('inf')
    #Here we try to calculate the time step to be taken on a particular iteration by ensuring the minimum possible value of delta_t
    for i in range(temperature_at_t.shape[0]):
        #The temperature of the t th time step is passed
        speed = speed_non_dimensional(temperature_at_t[i])
        del_t = courant_number*(deltax)/(speed + velocity_at_t[i])
        if(del_t<small):
            small = del_t
            
    return small


            
    

            

    
    