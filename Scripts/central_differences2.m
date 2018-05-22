function [ Q, t ] = central_differences2( gridpoints, initial_values, discrete_flux, timestepsize, finaltime)
%calculates the central difference scheme for the burger's equation
%equation 

%INPUT:
%gridpoints:        vector containing the grid points of the discretisation
%initial_values:    vector containing the initial values corresponding to
%the grid points
%discrete flux:     function handle for the discretisation of the flux
%timestepsize:      timestep size
%finaltime:         final time
%
%OUTPUT:
%Q:                 Matrix containing the numerical approximation of the
%solution, rows correspond to the time steps
%t:                 time vector

t = 0:timestepsize:finaltime;       %time vector
T = length(t);                      %number of time steps
k = length(gridpoints);             %number of gridpoints
Q = zeros(T, k+1);                  %initalize solution matrix

Q(1, :) = [initial_values(end) initial_values initial_values(1)];   %Initial values and periodic boundary conditions

for n = 1:T-1
    Q(n+1, 1) = Q(n, end-1);        %periodic boundary conditions
    Q(n+1, end) = Q(n, 2);          %periodic boundary conditions
    for j = 2:k
        a = Q(n, j);
        Q(n+1, j) = Q(n, j) - timestepsize/(gridpoints(j)-gridpoints(j-1))*( discrete_flux(Q(n,j), Q(n,j+1), a) - discrete_flux(Q(n,j-1), Q(n,j), a));    %calculation of the scheme
    end
end

end