import numpy
import numpy as np

def area(x):
    #Non dimensional
    # X lies between 0 and 3 
    return  1 + 2.2*(x-1.5)**2

def compute_diff(density,velocity,gamma,temperature,delta_x,a,a1,i):
        d1 = -density[i]*(velocity[i+1]-velocity[i])/(delta_x) - density[i]*velocity[i]*(np.log(a1)-np.log(a))/delta_x - velocity[i]*(density[i+1]-density[i])/delta_x
        d2 = -velocity[i]*(velocity[i+1]-velocity[i])/delta_x - (1/gamma)*((temperature[i+1]-temperature[i]))/delta_x + (temperature[i]/density[i])*(density[i+1]-density[i])/delta_x
        d3 = -velocity[i]*(temperature[i+1]-temperature[i])/delta_x - (gamma-1)*((temperature[i])/delta_x)*(velocity[i+1] - velocity[i] + velocity[i]*(np.log(a1)-np.log(a)))
        return d1,d2,d3
    
    

def diff_av_mccormack(Nx,density,temperature,velocity,delta_x,delta_t,gamma = 1.4):
    
    d_density,d_temperature,d_velocity = numpy.zeros_like(density),numpy.zeros_like(temperature),numpy.zeros_like(velocity)
    
    d_density_p,d_temperature_p,d_velocity_p = numpy.zeros_like(density),numpy.zeros_like(temperature),numpy.zeros_like(velocity)
    
    density_p,temperature_p,velocity_p = numpy.zeros_like(density),numpy.zeros_like(temperature),numpy.zeros_like(velocity) 
    
    d_density_av,d_temperature_av,d_velocity_av =numpy.zeros_like(density),numpy.zeros_like(temperature),numpy.zeros_like(velocity)
    
    
    numpy.seterr('raise')
    for i in range(Nx-2):
        a = area((i+1)*delta_x)
        a1 = area((i+2)*delta_x)
        d1,d2,d3 = compute_diff((density),(velocity),gamma,(temperature),delta_x,a,a1,i+1)
        d_density[i+1] = d1
        d_velocity[i+1] = d2
        d_temperature[i+1] = d3
        density_p[i+1] = density[i+1] + d1*delta_t
        velocity_p[i+1] = velocity[i+1] + d2*delta_t
        temperature_p[i+1] = temperature[i+1] + d3*delta_t
        
    for i in range(Nx-2):
        a = area((i+1)*delta_x)
        a1 = area((i+2)*delta_x)
        d1,d2,d3 = compute_diff(density_p,(velocity_p),gamma,(temperature_p),delta_x,a,a1,i+1)
        d_density_av[i+1] = 0.5*(d1 + d_density[i+1])
        d_velocity_av[i+1] = 0.5*(d2 + d_velocity[i+1])
        d_temperature_av[i+1] = 0.5*(d3 + d_temperature[i+1])

    
    return d_density_av,d_velocity_av,d_temperature_av