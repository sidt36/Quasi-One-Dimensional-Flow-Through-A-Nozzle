function[density,temperature,velocity] = boundary_condition(density,temperature,velocity)
    %The reservoir has only one floating variable , the inlet density
    %and temperature are fixed
    density(:,1) = 1;
    temperature(:,1) = 1;
    velocity(:,1) = 2.*velocity(:,2) - velocity(:,3);
    
    %All the flow variables are floating at the exit
    density(:,end) = 2.*density(:,end-1) - density(:,end-2);
    temperature(:,end) = 2.*temperature(:,end-1) - temperature(:,end-2);
    velocity(:,end) = 2.*velocity(:,end-1) - velocity(:,end-2);

    
end