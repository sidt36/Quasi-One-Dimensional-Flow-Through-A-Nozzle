
function [density,velocity,temperature,pressure,mach,mass_flow_rate_memory,pressure_throat,density_throat,velocity_throat,temperature_throat,Mach_throat,iter] = McCormack(Nx,x,delta_x,A,gamma,throat_index,courant_number)

%We initialize the data
[density,temperature,velocity] = initialize_flow_variables(Nx,1,delta_x);


change = inf;
tolerance = 1e-11;


  iter = 1;

% We march in time until our error metric is smaller than our tolerance 

  while iter<50000

    if change<tolerance
        break
    end
      
    density_old = density;
    velocity_old = velocity;
    temperature_old = temperature;
    mass_old = density_old.*A.*velocity_old;
    
     % Using the CFL criterion we predict the timestep to be taken
    [delta_t] = timestep(temperature,velocity,delta_x,courant_number);
    
    % Naive prediction
    
    for i = 2 : Nx-1
        
        diff1 = (density(i+1)-density(i))/delta_x;
        diff2 = (velocity(i+1)-velocity(i))/delta_x;
        diff3 = (log(A(i+1))-log(A(i)))/delta_x;
        diff4 = (temperature(i+1)-temperature(i))/delta_x;
        
        
        
        d_density(i) = (-density(i)*diff2) - ((density(i)*velocity(i))*diff3) - (velocity(i)*diff1);
        
        
        d_velocity(i) = (-velocity(i)*diff2) - ((1/gamma)*((diff4)+((temperature(i)/density(i))*diff1)));
        
        
        d_temperature(i) = (-velocity(i)*diff4) - ((gamma-1)*temperature(i))*(diff2 + (velocity(i)*diff3));    
        
        density(i) = density(i) + (d_density(i)*delta_t);
        velocity(i) = velocity(i) + (d_velocity(i)*delta_t);
        temperature(i) = temperature(i) + (d_temperature(i)*delta_t);
        
    end
    
    % We correct our naive prediction 
    for i = 2 : Nx-1
        
        diff1 = (density(i)-density(i-1))/delta_x;
        diff2 = (velocity(i)-velocity(i-1))/delta_x;
        diff3 = (temperature(i)-temperature(i-1))/delta_x;
        diff4 = (log(A(i))-log(A(i-1)))/delta_x;
        
        
        d_density_p(i) = (-density(i)*diff2) - ((density(i)*velocity(i))*diff4) - (velocity(i)*diff1);
        
        
        d_velocity_p(i) = (-velocity(i)*diff2) - ((1/gamma)*((diff3)+((temperature(i)/density(i))*diff1)));
        
        
        d_temperature_p(i) = (-velocity(i)*diff3) - ((gamma-1)*temperature(i))*(diff2 + (velocity(i)*diff4)); 
        
    end
    
    % We take the average
    
    d_density_avg = 0.5*(d_density+d_density_p);
    d_velocity_avg = 0.5*(d_velocity+d_velocity_p);
    d_temperature_avg = 0.5*(d_temperature+d_temperature_p);
    
    % Note the average arrays have size Nx-1 not Nx
    
 
    density = density_old + [d_density_avg*delta_t,0];
    velocity = velocity_old + [d_velocity_avg*delta_t,0];
    temperature = temperature_old  + [d_temperature_avg*delta_t,0];
    % We ensure the floating conditions are implemented
    
    [density,temperature,velocity] = boundary_condition(density,temperature,velocity);
    
    
    mach = velocity./sqrt(temperature);
    mass_flow = density.*A.*velocity;
    pressure = density.*temperature;
    
    % We moniter what happens in the throat
    
    density_throat(iter) = density(throat_index);
    temperature_throat(iter) = temperature(throat_index);
    Mach_throat(iter) = mach(throat_index);
    velocity_throat(iter) = velocity(throat_index);
    pressure_throat(iter) = pressure(throat_index);
    
    change = max(abs(mass_flow-mass_old));
    
    %We put our flow variables into memory  
    iter = iter+1;
 
    mass_flow_rate_memory(iter,:) = mass_flow;
    
    
       
end


end