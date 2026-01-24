
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This function is used to describe the parameter sampling of the pharmacokinetics of atezolizumab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S=PK_model_atezolizumab_parameter(sample,para)

%% Generating sampling matrix

S = zeros(sample,para);

S(:,1)=Beta_distribution(2.83,3.36,5,3,sample);       % V_1_A
S(:,2)=Beta_distribution(2.49,3.63,2,1,sample);       % V_2_A
S(:,3)=Beta_distribution(0.546,0.714,8,2,sample);     % Q_A
S(:,4)=Beta_distribution(0.2,0.274,9,1,sample);       % CI_A
S(:,5)=Beta_distribution(50,7000,7,2,sample);        % V_max_A
S(:,6)=Beta_distribution(1000,13000,9,2,sample);         % K_A

end

%% The mean and variance of the Gaussian distribution are calculated

% V_1_A = [2.57 2.99];
% a_V_1_A=(V_1_A(1)+V_1_A(2))/2;
% b_V_1_A=(V_1_A(2)-V_1_A(1))/6;
% 
% V_2_A = [2.29 2.57];
% a_V_2_A=(V_2_A(1)+V_2_A(2))/2;
% b_V_2_A=(V_2_A(2)-V_2_A(1))/6;
% 
% Q_A = [0.27 0.69];
% a_Q_A=(Q_A(1)+Q_A(2))/2;
% b_Q_A=(Q_A(2)-Q_A(1))/6;
% 
% CI_A = [0.15 0.30];
% a_CI_A=(CI_A(1)+CI_A(2))/2;
% b_CI_A=(CI_A(2)-CI_A(1))/6;
% 
% V_max = [120 130];
% a_V_max=(V_max(1)+V_max(2))/2;
% b_V_max=(V_max(2)-V_max(1))/6;
% 
% K_A = [50 100];
% a_K_A=(K_A(1)+K_A(2))/2;
% b_K_A=(K_A(2)-K_A(1))/6;

% S(:,1)=normrnd(a_V_1_A,b_V_1_A,sample,1);   % V_1_A
% S(:,2)=normrnd(a_V_2_A,b_V_2_A,sample,1);   % V_2_A
% S(:,3)=normrnd(a_Q_A,b_Q_A,sample,1);       % Q_A
% S(:,4)=normrnd(a_CI_A,b_CI_A,sample,1);     % CI_A
% S(:,5)=normrnd(a_V_max,b_V_max,sample,1);   % V_max
% S(:,6)=normrnd(a_K_A,b_K_A,sample,1);       % K_A

