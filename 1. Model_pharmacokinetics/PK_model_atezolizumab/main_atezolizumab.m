
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This script file is used to perform pharmacokinetic simulations of atezolizumab

% sample: the number of samples to perform the random simulation
% para: the number of parameters to perform the random simulation
% a: start time (therapy)
% b: end time (therapy)
% c: therapy period
% d: simulation time
% m: the number of therapy periods
% Drug_Ate: therapy dose of atezolizumab
% VP: Storage of numerical results on the pharmacokinetics of atezolizumab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sample = 1000;
para = 6;
a = 1;
d = 360;
c = 21;
Drug_Ate = 1200;

%% Parameter sampling

S = PK_model_atezolizumab_parameter(sample,para);
file_S='./VP/VP_PK_model_atezolizumab_parameter.mat';
save(file_S,'S');

%% Simulation calculation

[n,~]=size(S);

for i=1:n

    V_1_A = S(i,1);
    V_2_A = S(i,2);
    Q_A = S(i,3);
    CI_A = S(i,4);
    V_max = S(i,5);
    K_A = S(i,6);

    t_start = a;
    t_tau = c;
    t_end = t_start+t_tau;
    t_calibration = d;

    m=ceil(t_calibration/t_tau);

    Q_Ate_start = [Drug_Ate 0];
    Q_Ate = [Drug_Ate 0];

    tout = zeros(1,d);
    xout = zeros(d,2);

    for j = 1:m

        if j < m

            [t, x] = ode45(@(t,x) PK_model_atezolizumab(t,x,V_1_A,V_2_A,Q_A,CI_A,V_max,K_A),t_start:1:t_end,Q_Ate);

            nt = length(t);
            tout(t_start:1:t_end) = t_start:1:t_end;
            xout(t_start:1:t_end,:) = x(1:nt,:);
            
            t_start = t_end ;
            t_end = t_start + t_tau;
            Q_Ate = Q_Ate_start+x(end,:);

        else

            [t, x] = ode45(@(t,x) PK_model_atezolizumab(t,x,V_1_A,V_2_A,Q_A,CI_A,V_max,K_A),t_start:1:t_calibration,Q_Ate);

            nt = length(t);
            tout(t_start:1:t_calibration) = t_start:1:t_calibration;
            xout(t_start:1:t_calibration,:) = x(1:nt,:);

        end

    end

    xout(xout < 0.1 ) = 0.1; 

    file=['./VP/','VPsample',num2str(i),'.mat'];
    save(file,'xout');

end

%% Save results

VP_central=zeros(t_calibration,sample+1);
VP_peripheral=zeros(t_calibration,sample+1);

VP_central(:,1)=1:t_calibration;
VP_peripheral(:,1)=1:t_calibration;

for i=1:n

    file_sample=['./VP/','VPsample',num2str(i),'.mat'];
    A=cell2mat(struct2cell(load(file_sample)));
    VP_central(:,i+1)=A(:,1);
    VP_peripheral(:,i+1)=A(:,2);

end

file='./VP/VP_central.mat';
save(file,'VP_central');

file='./VP/VP_peripheral.mat';
save(file,'VP_peripheral');

