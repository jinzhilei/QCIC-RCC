
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This function is used to describe the parameter sampling of the pharmacokinetics of bevacizumab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S=PK_model_bevacizumab_parameter(sample,para)

%% Generating sampling matrix

S = zeros(sample,para);

S(:,1)=Beta_distribution(2.57,4.20,2,2,sample);       % V_1_B
S(:,2)=Beta_distribution(2.29,7.70,2,6,sample);       % V_2_B
S(:,3)=Beta_distribution(0.27,0.69,6,2,sample);       % Q_B
S(:,4)=Beta_distribution(0.15,0.30,1,9,sample);       % CI_B
S(:,5)=Beta_distribution(120,800,8,2,sample);         % V_max_B
S(:,6)=Beta_distribution(300,1300,4,5,sample);        % K_B

end

%% The mean and variance of the Gaussian distribution are calculated

% V_1_B = [2.57 2.99];
% a_V_1_B=(V_1_B(1)+V_1_B(2))/2;
% b_V_1_B=(V_1_B(2)-V_1_B(1))/6;
% 
% V_2_B = [2.29 2.57];
% a_V_2_B=(V_2_B(1)+V_2_B(2))/2;
% b_V_2_B=(V_2_B(2)-V_2_B(1))/6;
% 
% Q_B = [0.27 0.69];
% a_Q_B=(Q_B(1)+Q_B(2))/2;
% b_Q_B=(Q_B(2)-Q_B(1))/6;
% 
% CI_B = [0.15 0.30];
% a_CI_B=(CI_B(1)+CI_B(2))/2;
% b_CI_B=(CI_B(2)-CI_B(1))/6;
% 
% V_max = [120 130];
% a_V_max=(V_max(1)+V_max(2))/2;
% b_V_max=(V_max(2)-V_max(1))/6;
% 
% K_B = [50 100];
% a_K_B=(K_B(1)+K_B(2))/2;
% b_K_B=(K_B(2)-K_B(1))/6;

% S(:,1)=normrnd(a_V_1_B,b_V_1_B,sample,1);   % V_1_B
% S(:,2)=normrnd(a_V_2_B,b_V_2_B,sample,1);   % V_2_B
% S(:,3)=normrnd(a_Q_B,b_Q_B,sample,1);       % Q_B
% S(:,4)=normrnd(a_CI_B,b_CI_B,sample,1);     % CI_B
% S(:,5)=normrnd(a_V_max,b_V_max,sample,1);   % V_max
% S(:,6)=normrnd(a_K_B,b_K_B,sample,1);       % K_B

