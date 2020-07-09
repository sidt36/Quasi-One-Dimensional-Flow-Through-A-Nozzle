function[density,temperature,velocity] = initialize_flow_variables(spacesteps,timesteps,delta_x)
    %Create Placeholder 
    density = zeros(timesteps,spacesteps);
    temperature = zeros(timesteps,spacesteps);
    velocity = zeros(timesteps,spacesteps);
    
    %Create intial values which help with better convergence
    
    for i = 1:spacesteps
        x = (i-1)*delta_x;
        density(1,i) = 1-(0.3146*x); 
        temperature(1,i) = 1-(0.2314*x);
        velocity(1,i) = (0.1+(1.09*x)).*(temperature(1,i)).^0.5;
        
    end
    
end