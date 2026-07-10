function ROC()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  This script is used to calculate the ROC value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sample = 1664;
r = zeros(sample,1);
ratio = zeros(sample,1);
Imm = zeros(sample,32);
day = 84;
for i =1:sample
    file = ['../RECISE/4. AB/VP_AB/VPsample_',num2str(i),'.mat'];
    Tumor_0 = cell2mat(struct2cell(load(file)));
    cell_0 = ((3 .* Tumor_0(1,25) .* 2.572e-9 .* 0.8) ./ (0.37 * 4 * 3.14)).^(1/3);
    cell = ((3 .* (Tumor_0(day,25)+Tumor_0(day,26)+Tumor_0(day,27)+Tumor_0(day,28)) .* 2.572e-9 .* 0.8) ./ (0.37 * 4 * 3.14)).^(1/3);
    ratio(i) = (cell - cell_0) ./ cell_0.*100;
    Imm(i,:) = Tumor_0(day,:);
end

for i = 1:sample
        if ratio(i) >= 20
            r(i) = 0;
        else 
            r(i) = 1;
        end
end
save('./VP_immune/response.mat',"r");
save('./VP_immune/immune.mat','Imm')

end