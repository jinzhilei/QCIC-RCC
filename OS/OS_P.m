function OS_P()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This script is used to calculate the OS of patients under different indicators

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sample = 1664;
day = 84;
T = zeros(7,sample);
for i =1:sample
    file = ['../4. AB/VP_AB/VPsample_',num2str(i),'.mat'];
    Tumor_0 = cell2mat(struct2cell(load(file)));
    T(1,i) = Tumor_0(day,12);
    T(2,i) = Tumor_0(day,15);
    T(3,i) = Tumor_0(day,21);
    T(4,i) = Tumor_0(day,24);
end
T(5,:) = T(2,:)./T(3,:);
T(6,:) = T(3,:)./T(4,:);
T(7,:) = T(1,:)./T(3,:);

T_m = median(T,2);
Num = cell(8, 1);
rows_all = 1:sample;
L_1 = find(T(3,:)>T_m(3));
rows_other_1 = setdiff(rows_all, L_1);
L_2 = find(T(5,:)>T_m(5));
rows_other_2 = setdiff(rows_all, L_2);
L_3 = find(T(6,:)>T_m(6));
rows_other_3 = setdiff(rows_all, L_3);
L_4 = find(T(7,:)>T_m(7));
rows_other_4 = setdiff(rows_all, L_1);
Num{1, 1} = L_1;
Num{2, 1} = rows_other_1;
Num{3, 1} = L_2;
Num{4, 1} = rows_other_2;
Num{5, 1} = L_3;
Num{6, 1} = rows_other_3;
Num{7, 1} = L_4;
Num{8, 1} = rows_other_4;

file = './data/death_all.mat';
death_all = cell2mat(struct2cell(load(file)));
for a =1:8
    death = death_all(Num{a,1},:);
    n = length(Num{a,1});
    OS(1)=n;
    for j=2:1:360
        OS(j)=OS(j-1)-length(find(death==j))*1;
        if OS(j)<0
            OS(j)=0;
        end
    end
    Q_ave=OS./n*100;

    L = randi([1,n],100,100);
    for i =1:100
        Death=sort(death(L(:,i)));

        OS(1)=100;
        for j=2:1:360
            OS(j)=OS(j-1)-length(find(Death==j))*1;
            if OS(j)<0
                OS(j)=0;
            end
        end
        Q(:,i)=OS;
    end
    file = ['./data/OS_ave',num2str(a),'.mat'];
    save(file,"Q_ave");
    file = ['./data/OS',num2str(a),'.mat'];
    save(file,"Q");
end


end