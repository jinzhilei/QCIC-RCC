function solve()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Calculate the ROC value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sample = 10000;
day = 360;
jianc =  43;
tumor_NoDrug = zeros(day,sample);
tumor_A = zeros(day,sample);
tumor_AB = zeros(day,sample);
ratio_NoDrug = zeros(day,sample);
ratio_A = zeros(day,sample);
ratio_AB = zeros(day,sample);
A = zeros(3,sample);
B = zeros(3,sample);
file_R_NoDrug='../5.VP/VP_NoDrug/VP_NoDrug_bar/VP_tumor_number_roc.mat';
r_NoDrug=cell2mat(struct2cell(load(file_R_NoDrug)));
tumor_0_NoDrug = ((r_NoDrug.*2.572e-9.*0.8*3)./(0.37*4*3.14)).^(1/3);

file_R_A='../5.VP/VP_A/VP_A_bar/VP_tumor_initialzation_roc.mat';
r_A=cell2mat(struct2cell(load(file_R_A)));
tumor_0_A = ((r_A.*2.572e-9.*0.8*3)./(0.37*4*3.14)).^(1/3);

file_R_AB='../5.VP/VP_AB/VP_AB_bar/VP_tumor_initialzation_roc.mat';
r_AB=cell2mat(struct2cell(load(file_R_AB)));
tumor_0_AB = ((r_AB.*2.572e-9.*0.8*3)./(0.37*4*3.14)).^(1/3);

%% Calculation of tumor evaluation indicators

for i =1:sample
    file_1=['../5.VP/VP_NoDrug/VP_NoDrug_bar_number/','VPsample',num2str(i),'.mat'];
    file_2=['../5.VP/VP_A/VP_A_bar_number/','VPsample_A',num2str(i),'.mat'];
    file_3=['../5.VP/VP_AB/VP_AB_bar_number/','VPsample_AB',num2str(i),'.mat'];
    cell_NoDrug=cell2mat(struct2cell(load(file_1)));
    cell_A=cell2mat(struct2cell(load(file_2)));
    cell_AB=cell2mat(struct2cell(load(file_3)));
    tumor_NoDrug(:,i) = cell_NoDrug(:,25)+cell_NoDrug(:,26)+cell_NoDrug(:,27)+cell_NoDrug(:,28);
    tumor_NoDrug(:,i) = ((3*tumor_NoDrug(:,i).*2.572e-9.*0.8)./(0.37*4*3.14)).^(1/3);
    ratio_NoDrug(:,i) = (tumor_NoDrug(:,i)- tumor_0_NoDrug(i))./tumor_0_NoDrug(i)*100;
    tumor_A(:,i) = cell_A(:,25)+cell_A(:,26)+cell_A(:,27)+cell_A(:,28);
    tumor_A(:,i) = ((3*tumor_A(:,i).*2.572e-9.*0.8)./(0.37*4*3.14)).^(1/3);
    ratio_A(:,i) = (tumor_A(:,i)- tumor_0_A(i))./tumor_0_A(i)*100;
    tumor_AB(:,i) = cell_AB(:,25)+cell_AB(:,26)+cell_AB(:,27)+cell_AB(:,28);
    tumor_AB(:,i) = ((3*tumor_AB(:,i).*2.572e-9.*0.8)./(0.37*4*3.14)).^(1/3);
    ratio_AB(:,i) = (tumor_AB(:,i)- tumor_0_AB(i))./tumor_0_AB(i)*100;
end

%% Plot

for i = 1:sample
    if ratio_NoDrug(day,i) >= 20
        A(1,i) = 1;
    elseif -30 <= ratio_NoDrug(day,i) && ratio_NoDrug(day,i) < 20
        A(1,i) = 2;
    elseif  -80 <= ratio_NoDrug(day,i) &&ratio_NoDrug(day,i) <-30
       A(1,i) = 3;
    elseif  ratio_NoDrug(day,i) <-80
       A(1,i) = 4;
    end    

    if ratio_NoDrug(jianc:day,i) >= -30
        B(1,i) = 0;
    elseif  ratio_NoDrug(jianc:day,i) <-30
        B(1,i) = 1;
    end    

end

for i = 1:sample
    if ratio_A(day,i) >= 20
        A(2,i) = 1;
    elseif -30 <= ratio_A(day,i) && ratio_A(day,i) < 20
        A(2,i) = 2;
    elseif  -80 <= ratio_A(day,i) &&ratio_A(day,i) <-30
       A(2,i) = 3;
    elseif  ratio_A(day,i) <-80
       A(2,i) = 4;
    end    

    if ratio_NoDrug(jianc:day,i) >= -30
        B(2,i) = 0;
    elseif  ratio_NoDrug(jianc:day,i) <-30
        B(2,i) = 1;
    end    

end

for i = 1:sample
    if ratio_AB(day,i) >= 20
        A(3,i) = 1;
    elseif -30 <= ratio_AB(day,i) && ratio_AB(day,i) < 20
        A(3,i) = 2;
    elseif  -80 <= ratio_AB(day,i) &&ratio_AB(day,i) <-30
       A(3,i) = 3;
    elseif  ratio_AB(day,i) <-80
       A(3,i) = 4;
    end    

    if ratio_AB(jianc:day,i) >= -30
        B(3,i) = 0;
    elseif  ratio_AB(jianc:day,i) <-30
        B(3,i) = 1;
    end    

end

file_sort='./data/immune_sort.mat';
save(file_sort,'A')

file_response='./data/immune_response.mat';
save(file_response,'B')

end