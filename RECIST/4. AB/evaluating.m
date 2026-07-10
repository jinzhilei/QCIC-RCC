function [CR, PR, SD, PD] = evaluating(sample, you)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This script file is used to perform the assessment virtual patient

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    file = './VP_immune/Tumor_0';
    Tumor_0 = cell2mat(struct2cell(load(file)));
    cell_0 = ((3 .* Tumor_0 .* 2.572e-9 .* 0.8) ./ (0.37 * 4 * 3.14)).^(1/3);
    cell = ((3 .* (you(:,25)+you(:,26)+you(:,27)+you(:,28)) .* 2.572e-9 .* 0.8) ./ (0.37 * 4 * 3.14)).^(1/3);
    ratio = (cell - cell_0) ./ cell_0.*100;
    
    A = zeros(1, 4);
    for i = 1:sample
        if ratio(i) >= 20
            A(1) = A(1) + 1;
        elseif -30 <= ratio(i) && ratio(i) < 20
            A(2) = A(2) + 1;
        elseif -80 <= ratio(i) && ratio(i) < -30
            A(3) = A(3) + 1;
        elseif ratio(i) < -80
            A(4) = A(4) + 1;
        end
    end
    
    PD = A(1) / sample * 100;
    SD = A(2) / sample * 100;
    PR = A(3) / sample * 100;
    CR = A(4) / sample * 100;
end