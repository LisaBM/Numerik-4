discrete_flux = @(left_value, right_value, a)a/2*(left_value^2 + right_value^2);
timestepsize = 1/100;
finaltime = timestepsize * 50;
k = 100;
gridpoints = 0:1/k:2*pi;
midpoints = 1/2*(gridpoints(1:end-1) + gridpoints(2:end));
%initial values
q1 = @(x)sin(x);
initial_values1 = q1(midpoints);

%calculation of the difference schemes
[ Q1, t ] = central_differences2( gridpoints, initial_values1, discrete_flux, timestepsize, finaltime );
[ Q2, ~ ] = upwind_approximation( gridpoints, initial_values1, discrete_flux, timestepsize, finaltime );

figure
for i=1:size(Q1, 1)
    plot(midpoints, Q1(i, 2:end-1), midpoints, Q2(i, 2:end-1))
    ylim([-2 2])
    legend('central differences', 'upwind')
    pause(0.1)
end
