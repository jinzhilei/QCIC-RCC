function  S = LHS(n_samples, n_dims, bounds)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% This script file is used to extract parameters for Latin hypercube extraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    rng(123);

    S = zeros(n_samples, n_dims);
    for i = 1:n_dims
        cut = linspace(0, 1, n_samples + 1);
        for j = 1:n_samples
            S(j, i) = cut(j) + rand * (cut(j+1) - cut(j));
        end
        S(:, i) = S(randperm(n_samples), i);
    end
    if ~isempty(bounds)
        for i = 1:n_dims
            low = bounds(i, 1);
            high = bounds(i, 2);
            S(:, i) = low + S(:, i) * (high - low);
        end
    end
end