function OS()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This script is used to calculate the OS
% sample: Number of virtual patients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sample = 1664;
OS = zeros(1,360);
Q = zeros(360,100);
param_bounds =[0.1,1;0.05,0.5];
dim = length(param_bounds);

for i =1:sample
file = ['../RECIST/4. AB/VP_AB/VPsample_',num2str(i),'.mat'];
Tumor = cell2mat(struct2cell(load(file)));
cell(:,4*i-3:4*i) = Tumor(:,25:28);
end

file = '../4. AB/VP_immune/Tumor_cell.mat';
G = cell2mat(struct2cell(load(file)));
X_1 = load('data/1.Rini.2019.dat');

params1 = [0; cumsum(sqrt(diff(X_1(:,1)).^2 + diff(X_1(:,2)).^2))];
params1 = params1 / params1(end);
dense_n = 1000; 
params_dense = linspace(0, 1, dense_n)';
x1_dense = spline(params1, X_1(:,1), params_dense);
y1_dense = spline(params1, X_1(:,2), params_dense);

target_dense = [x1_dense, y1_dense];

plot(x1_dense, y1_dense);
hold on
plot(X_1(:,1),X_1(:,2))

X = zeros(360,sample);
  for i =1:1664
   X(:,i) = 0.1*cell(:,4*i-3)./G(i,1)+0.2*cell(:,4*i-2)./G(i,4)+0.2*cell(:,4*i-1)./G(i,9)+0.5*cell(:,4*i)./G(i,10);
  end

pool = gcp('nocreate');
if ~isempty(pool)
    delete(pool);
end
parpool(20);  
pool = gcp;

best_params = PSO(X,target_dense,param_bounds,dim,pool);
save('data/best.mat',"best_params");

P_death=1./(1+exp(-(X-best_params(1))./best_params(2)));
P = rand(size(P_death));
deathEvent = P < P_death;
[row, col] = find(deathEvent);
Death_time = accumarray(col, row, [size(P_death,2), 1], @min);
Death_all=sort(Death_time);
save('./data/death_all.mat',"Death_time")
OS(1)=1664;
    for j=2:1:360
        OS(j)=OS(j-1)-length(find(Death_all==j))*1;
        if OS(j)<0
           OS(j)=0;
        end
    end
    Q_ave=OS./1664*100;

L = randi([1,sample],100,100);
for i =1:100
P_death=1./(1+exp(-(X(:,L(:,i))-best_params(1))./best_params(2)));
P = rand(size(P_death));
deathEvent = P < P_death;
[row, col] = find(deathEvent);
Death_time = accumarray(col, row, [size(P_death,2), 1], @min);
Death=sort(Death_time);

OS(1)=100;
    for j=2:1:360
        OS(j)=OS(j-1)-length(find(Death==j))*1;
        if OS(j)<0
           OS(j)=0;
        end
    end
    Q(:,i)=OS;
end
save('./data/OS_ave.mat',"Q_ave");
save('./data/OS.mat',"Q");



for i =1:sample
file = strcat('./treat_AB/treat2/(2,2)_VPsample_AB',num2str(i),'.mat');
Tumor = cell2mat(struct2cell(load(file)));
cell(:,4*i-3:4*i) = Tumor(:,25:28);
end
X = zeros(360,sample);
  for i =1:1664
   X(:,i) = 0.1*cell(:,4*i-3)./G(i,1)+0.2*cell(:,4*i-2)./G(i,4)+0.2*cell(:,4*i-1)./G(i,9)+0.5*cell(:,4*i)./G(i,10);
  end


P_death=1./(1+exp(-(X-best_params(1))./best_params(2)));
P = rand(size(P_death));
deathEvent = P < P_death;
[row, col] = find(deathEvent);
Death_time = accumarray(col, row, [size(P_death,2), 1], @min);
Death=sort(Death_time);

OS(1)=1664;
    for j=2:1:360
        OS(j)=OS(j-1)-length(find(Death==j))*1;
        if OS(j)<0
           OS(j)=0;
        end
    end
    Q_ave=OS./1664*100;
L = randi([1,sample],100,100);
for i =1:100
P_death=1./(1+exp(-(X(:,L(:,i))-best_params(1))./best_params(2)));
P = rand(size(P_death));
deathEvent = P < P_death;
[row, col] = find(deathEvent);
Death_time = accumarray(col, row, [size(P_death,2), 1], @min);
Death=sort(Death_time);

OS(1)=100;
    for j=2:1:360
        OS(j)=OS(j-1)-length(find(Death==j))*1;
        if OS(j)<0
           OS(j)=0;
        end
    end
    Q(:,i)=OS;
end
save('./data/OS_ave_2.mat',"Q_ave");
save('./data/OS_2.mat',"Q");


for i =1:sample
file = strcat('./treat_AB/treat3/(2,2)_VPsample_AB',num2str(i),'.mat');
Tumor = cell2mat(struct2cell(load(file)));
cell(:,4*i-3:4*i) = Tumor(:,25:28);
end
X = zeros(360,sample);
  for i =1:1664
   X(:,i) = 0.1*cell(:,4*i-3)./G(i,1)+0.2*cell(:,4*i-2)./G(i,4)+0.2*cell(:,4*i-1)./G(i,9)+0.5*cell(:,4*i)./G(i,10);
  end


P_death=1./(1+exp(-(X-best_params(1))./best_params(2)));
P = rand(size(P_death));
deathEvent = P < P_death;
[row, col] = find(deathEvent);
Death_time = accumarray(col, row, [size(P_death,2), 1], @min);
Death=sort(Death_time);

OS(1)=1664;
    for j=2:1:360
        OS(j)=OS(j-1)-length(find(Death==j))*1;
        if OS(j)<0
           OS(j)=0;
        end
    end
    Q_ave=OS./1664*100;
L = randi([1,sample],100,100);
for i =1:100
P_death=1./(1+exp(-(X(:,L(:,i))-best_params(1))./best_params(2)));
P = rand(size(P_death));
deathEvent = P < P_death;
[row, col] = find(deathEvent);
Death_time = accumarray(col, row, [size(P_death,2), 1], @min);
Death=sort(Death_time);

OS(1)=100;
    for j=2:1:360
        OS(j)=OS(j-1)-length(find(Death==j))*1;
        if OS(j)<0
           OS(j)=0;
        end
    end
    Q(:,i)=OS;
end
save('./data/OS_ave_3.mat',"Q_ave");
save('./data/OS_3.mat',"Q");


end