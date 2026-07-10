function counts_recurrence_rate()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This script file is used to perform the cell dynamics simulation of different drug treatments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sample =100;
drug_number = 4;
A = zeros(drug_number,drug_number);
tumor_0 = zeros(sample,1);
for p = 1:3
    disp(p)
    for m =1:drug_number
        for n = 1:drug_number
            for i =1:sample
                file=strcat('treat_AB/treat',num2str(p),'/(',num2str(m),',',num2str(n),')','_VPsample_AB',num2str(i),'.mat');
                cell=cell2mat(struct2cell(load(file)));
                tumor_0(i) = ((3*cell(1,25).*2.572e-9.*0.8)./(0.37*4*3.14)).^(1/3);
                tumor(:,i) = cell(:,25)+cell(:,26)+cell(:,27)+cell(:,28);
                tumor(:,i) = ((3*tumor(:,i).*2.572e-9.*0.8)./(0.37*4*3.14)).^(1/3);
                ratio(:,i) = (tumor(:,i)- tumor_0(i))/tumor_0(i);
                if ratio(84,i) >= 0.2                    
                    A(m,n) = A(m,n)+1;
                end

            end
        end
    end
    A = A./sample;
    file=['Recurrence_rate/treat',num2str(p),'_Recurrence_rate_',num2str(84),'.mat'];
    save(file,'A');
end
end

