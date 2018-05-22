A = [-1, 1];
L = [0.5, 1, 3/2];
K = [10, 100]; 

q1 = @(x)sin(x);
q2 = @(x) (x>=1/3 & x<=2/3);
        
for index_A = 1:2 
    for index_L = 1:3
        for index_K = 1:2
            %set  parameters
            lambda = L(index_L);
            k = K(index_K);
            title_text = strcat('a=', num2str(a), ', k=', num2str(k), ', lambda=', num2str(lambda));
            discrete_flux1 = @(left_value, right_value, a)a/2*(left_value + right_value);
            discrete_flux2 = @(left_value, right_value, a)a*left_value;
            timestepsize = lambda*1/k;
            finaltime = timestepsize * 50;
            gridpoints = 0:1/k:2*pi;
            midpoints = 1/2*(gridpoints(1:end-1) + gridpoints(2:end));
            a = A(index_A);
            %initial values
            initial_values1 = q1(midpoints);
            initial_values2 = q2(midpoints);
            %calculation of the difference schemes
            [Q1, t] = central_differences(gridpoints, initial_values1, discrete_flux1, timestepsize, finaltime, a);
            [Q2, ~] = central_differences(gridpoints, initial_values2, discrete_flux1, timestepsize, finaltime, a);
            [Q3, ~] = central_differences(gridpoints, initial_values1, discrete_flux2, timestepsize, finaltime, a);
            [Q4, ~] = central_differences(gridpoints, initial_values2, discrete_flux2, timestepsize, finaltime, a);
            true_solution1 = @(x, t)q1(x-a*t);
            true_solution2 = @(x, t)q2(x-a*t);
            %plots
            figure('units','normalized','outerposition',[0 0 1 1])
            for i = 1:size(Q1, 1)
                subplot(1,2,1)
                plot(midpoints, true_solution1(midpoints, t(i)), midpoints, Q1(i, 2:end-1), midpoints, Q3(i, 2:end-1))
                legend('true solution', 'numerical approximation w/ flux 1', 'numerical approximation w/ flux 2')
                title(strcat('initial values 1 ,', title_text))
                ylim([-2 2])

                subplot(1,2,2)
                plot(midpoints, true_solution2(midpoints, t(i)), midpoints, Q2(i, 2:end-1), midpoints, Q4(i, 2:end-1))
                legend('true solution', 'numerical approximation w/ flux 1', 'numerical approximation w/ flux 2')
                title(strcat('initial values 2 ,', title_text))
                ylim([-1 2])

                pause(0.1)               
            end
        end
    end
end