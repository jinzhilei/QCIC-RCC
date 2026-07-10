function B = e(sample)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This script file is used to perform the assessment virtual patient

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ratio = zeros(sample,1);
    day = 84;
for i =1:sample
    file = ['./VP_A/VPsample_',num2str(i),'.mat'];
    Tumor_0 = cell2mat(struct2cell(load(file)));
    cell_0 = ((3 .* Tumor_0(1,25) .* 2.572e-9 .* 0.8) ./ (0.37 * 4 * 3.14)).^(1/3);
    cell = ((3 .* (Tumor_0(day,25)+Tumor_0(day,26)) .* 2.572e-9 .* 0.8) ./ (0.37 * 4 * 3.14)).^(1/3);
    ratio(i) = (cell - cell_0) ./ cell_0.*100;
end
r = randi([1 sample], 100, 10);
    A = zeros(10, 4);
    B = zeros(1,4);
    for j = 1:10
    for i = 1:100
        if ratio(r(i,j)) >= 20
            A(j,1) = A(j,1) + 1;
        elseif -30 <= ratio(r(i,j)) && ratio(r(i,j)) < 20
            A(j,2) = A(j,2) + 1;
        elseif -80 <= ratio(r(i,j)) && ratio(r(i,j)) < -30
            A(j,3) = A(j,3) + 1;
        elseif ratio(r(i,j)) < -80
            A(j,4) = A(j,4) + 1;
        end
    end
    end
    for i = 1:sample
        if (ratio(i)) >= 20
            B(1) = B(1) + 1;
        elseif -30 <= (ratio(i)) && (ratio(i)) < 20
            B(2) = B(2) + 1;
        elseif -80 <= (ratio(i)) && (ratio(i)) < -30
            B(3) = B(3) + 1;
        elseif (ratio(i)) < -80
            B(4) = B(4) + 1;
        end
    end
    B = B./sample*100;
    save('./VP_immune/CPSP.mat',"A")    
    save('./VP_immune/CPSP_Ave.mat',"B")

end