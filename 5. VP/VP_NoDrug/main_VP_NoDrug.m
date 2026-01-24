
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This script file is used to perform cell dynamics simulations of no drug

% t_start: start time (therapy)
% t_calibration: simulation time
% z: simulate the number of virtual patients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_start = 1;
t_calibration = 360;
z = 100;

%% Parameter sampling

file='./VP_patient/VP_model_parameter.mat';
A_immune_par=cell2mat(struct2cell(load(file)));

file = './VP_patient/Number_VP_patients.mat';
O_order=cell2mat(struct2cell(load(file)));
L  = length(O_order);

%% Calculation of tumor parameters

S = tumor_parameter(z);
file_S='./VP_patient/VP_tumor_X1_parameter.mat';
save(file_S,'S');

B_immune_par = A_immune_par(randi([1 L], z, 1),:);

%% Simulation calculation

global par
r = 3.2279e10+rand(z,1) *(2.0675e11-3.2279e10);
file_R='./VP_NoDrug/VP_tumor_number.mat';
save(file_R,'r');
for i=1:z
    disp(i)

    theta_TN4 = B_immune_par(i,1);
    theta_TN8 = B_immune_par(i,2);
    theta_D0 = B_immune_par(i,3);
    theta_M = B_immune_par(i,4);
    kappa_Th = B_immune_par(i,5);
    kappa_Tr = B_immune_par(i,6);
    kappa_Tc = B_immune_par(i,7);
    beta_Th = B_immune_par(i,8);
    beta_Tr = B_immune_par(i,9);
    beta_Tc = B_immune_par(i,10);

    G_X1 = S(i,1);
    beta_X1 = S(i,2);
    d_X1 = S(i,3);

    x0=initialzation_VP_NoDrug(r(i));
    [t,y] = ode45(@(t,x) QCIC_VP_NoDrug(t,x,theta_TN4,theta_TN8,theta_D0,theta_M,kappa_Th,kappa_Tr,kappa_Tc,beta_Th,beta_Tr,beta_Tc,G_X1,beta_X1,d_X1),t_start:1:t_calibration, x0);
    nt = length(t);
    tout(t_start:1:t_calibration) = t_start:1:t_calibration;
    you(t_start:1:t_calibration,:) = y(1:nt,:);

    file=['./VP_NoDrug/','VPsample_',num2str(i),'.mat'];
    save(file,'you');
end
evaluating(z)
% picture_bar(z)


%% This function is used to describe the parameter sampling of VP

function S = tumor_parameter(z)
S = zeros(z,3);

S(:,1) = Beta_distribution(5e12,5e13,9,3,z);       % G_X1         %1.5e9
S(:,2) = Beta_distribution(0.1,0.2,90,90,z);      % \beta_X1     %0.15
S(:,3) = Beta_distribution(0.03,0.15,50,20,z);      % \d_X1       %0.04
% S(:,2) = Beta_distribution(0.1,0.2,65,50,z);      % \beta_X1     %0.15
% S(:,3) = Beta_distribution(0.02,0.2,40,50,z);      % \d_X1       %0.04



end

 
%% This function is used to generate random numbers that follow a Beta distribution

function y=Beta_distribution(a,b,c,d,e)

y=a+(b-a)*betarnd(c,d,[1,e]);

end