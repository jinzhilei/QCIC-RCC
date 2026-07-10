function [best_params] = PSO(Par, dim, param_bounds, t_start, t_calibration,pool)
%% This script is used to perform PSO calculations

% n_particles : Number of particles
% n_iterations: Number of iterations
% w: Inertia weight
% w_damp: Inertia weight attenuation coefficient
% c1: Individual learning factor
%c2: Social learning factor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n_particles = 20;      
    n_iterations = 10;           
    
    w = 0.4;       
    w_damp = 0.99;          
    c1 = 0.5;          
    c2 = 2.5;             
    
    
    %% Initialize particle swarm
    positions = zeros(n_particles, dim);
    for i = 1:dim
        positions(:,i) = param_bounds(i,1) + (param_bounds(i,2) - param_bounds(i,1)) * rand(n_particles, 1);
    end
    
    velocities = 0.1 * (param_bounds(:,2)' - param_bounds(:,1)') .* randn(n_particles, dim);
    
    pbest_errors = zeros(1, n_particles);
    futures = parallel.Future.empty();
    
    for i = 1:n_particles
        futures(i) = parfeval(pool, @calculate_error_parallel, 1, positions(i,:), Par, t_start, t_calibration);
    end
    
    for i = 1:n_particles
        [idx, error] = fetchNext(futures);
        pbest_errors(idx) = error;
    end
    
    pbest_positions = positions;
    

    [gbest_error, gbest_index] = min(pbest_errors);
    gbest_position = pbest_positions(gbest_index, :);
   
    history = zeros(n_iterations, 2);
    
    %% Main cycle
    for iter = 1:n_iterations

        new_positions = positions;
        new_velocities = velocities;
        
        for i = 1:n_particles
            r1 = rand(1, dim);
            r2 = rand(1, dim);
            new_velocities(i,:) = w * velocities(i,:) ...
                + c1 * r1 .* (pbest_positions(i,:) - positions(i,:)) ...
                + c2 * r2 .* (gbest_position - positions(i,:));
            
            new_positions(i,:) = positions(i,:) + new_velocities(i,:);
            
            for j = 1:dim
                new_positions(i,j) = max(new_positions(i,j), param_bounds(j,1));
                new_positions(i,j) = min(new_positions(i,j), param_bounds(j,2));
                
                if new_positions(i,j) == param_bounds(j,1) || new_positions(i,j) == param_bounds(j,2)
                    new_velocities(i,j) = 0;
                end
            end
        end

        futures = parallel.Future.empty();
        for i = 1:n_particles
            futures(i) = parfeval(pool, @calculate_error_parallel, 1, ...
                new_positions(i,:), Par, t_start, t_calibration);
        end
        

        for i = 1:n_particles
            [idx, current_error] = fetchNext(futures);
            
            if current_error < pbest_errors(idx)
                pbest_errors(idx) = current_error;
                pbest_positions(idx,:) = new_positions(idx,:);
            end
            
            if current_error < gbest_error
                gbest_error = current_error;
                gbest_position = new_positions(idx,:);
            end
        end
        
        positions = new_positions;
        velocities = new_velocities;
        
        w = w * w_damp;
        
        history(iter, 1) = gbest_error;
        history(iter, 2) = w;
    end
    file = './VP_immune/history';
    save(file,"history");
    
    best_params = gbest_position;
    
    %% Draw convergence curve
    figure('Position', [100, 100, 1200, 500]);
    plot(1:n_iterations, history(:,1), 'b-', 'LineWidth', 2);
    xlabel('Iteration count');
    ylabel('Loss function value');
    title('PSO convergence curve');
    grid on;
    set(gca, 'YScale', 'log');
    set(gca,'FontName','Arial')
end

function error = calculate_error_parallel(p, Par, t_start, t_calibration)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This script is used to calculate the error function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    m = size(Par, 1);
    S = tumor_parameter(m, p);
    you = zeros(m, 28);

    pool = gcp('nocreate');
    if isempty(pool)
        for i = 1:m
            y = solve_one_patient(i,Par, S, t_start, t_calibration);
            you(i,:) = y(end,:);
        end
    else
        futures = parallel.Future.empty();
        for i = 1:m
            futures(i) = parfeval(pool, @solve_one_patient, 1, i, Par, S, t_start, t_calibration);
        end
        
        for i = 1:m
            [~, result] = fetchNext(futures);
            you(i,:) = result(end,:);
        end
    end
    
    [CR, PR, SD, PD] = evaluating(m, you);
    v = [1 1 1 1];
  data_v = abs([CR-0,PR-8.85,SD-57.1,PD-34.05]');
  error = v * data_v;

end

