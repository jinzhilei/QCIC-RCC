function evaluating()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This module is used to describe the tumor evaluation metrics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameter sampling
N =10;
sample = 100;
jianc = 84;
CR = zeros(N,3);
PR = zeros(N,3);
SD = zeros(N,3);
PD = zeros(N,3);

Num = randi(1000,10,sample);

for p =1:3
    file=strcat('Recurrence_rate/treat',num2str(p),'_ratio(',num2str(2),',',num2str(2),').mat');
    R=cell2mat(struct2cell(load(file)));
    for j = 1:N
        for i =1:sample
            if R(jianc,Num(j,i)) >= 0.2
                PD(j,p) = PD(j,p)+1;
            elseif -0.3 <= R(jianc,Num(j,i)) && R(jianc,Num(j,i)) < 0.2
                SD(j,p) = SD(j,p)+1;
            elseif  -0.8 <= R(jianc,Num(j,i)) && R(jianc,Num(j,i)) <-0.3
                PR(j,p) = PR(j,p)+1;
            elseif  R(jianc,Num(j,i)) <-0.8
                CR(j,p) = CR(j,p)+1;
            end    
        end
    end
end
PD = PD./sample;SD = SD./sample;PR = PR./sample;CR = CR./sample;
file='Recurrence_rate/PD.mat';
save(file,'PD');
file='Recurrence_rate/SD.mat';
save(file,'SD');
file='Recurrence_rate/PR.mat';
save(file,'PR');
file='Recurrence_rate/SR.mat';
save(file,'CR');

end