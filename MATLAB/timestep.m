%We try to calculate the timestep by which each iteration should march by 
function[delta_t] = timestep(temperature_at_t,velocity_at_t,delta_x,C)
    
    speed = sqrt(temperature_at_t);
    
    delta_T = C*(delta_x)./(velocity_at_t + speed);
    
    delta_t = min(delta_T);
    
    
    
end