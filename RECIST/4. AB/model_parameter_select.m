function samples = model_parameter_select(sample, pare)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This script file is used to perform immune parameter sampling

% sample: Number of parameter samples
% pare: Number of parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    file_index = '../../VP_Immunogenomic/VP_immune/PDF_immune.mat';
    PDF = cell2mat(struct2cell(load(file_index)));
    

    samples = zeros(sample, pare);
    
    for i = 1:pare
        x = PDF(2*i-1, :);
        f = PDF(2*i,   :);
        
        cdf = cumtrapz(x, f);
        cdf = cdf / cdf(end);   
        
        r = rand(sample, 1);
        samples(:, i) = interp1(cdf, x, r, 'linear', 'extrap');
    end
end