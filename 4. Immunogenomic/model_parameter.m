%% This function is used to describe the parameter sampling of VPs

function S=model_parameter(sample,pare)

S = zeros(sample,pare);

S(:,1)=Beta_distribution(5e6,8e9,1,4,sample);         % \theta_T_N4
S(:,2)=Beta_distribution(5e6,8e9,8,4,sample);         % \theta_T_N8
S(:,3)=Beta_distribution(1e7,1e9,1,9,sample);         % \theta_D_0
S(:,4)=Beta_distribution(1e7,8e8,.1,9,sample);         % \theta_M 
S(:,5)=Beta_distribution(0.1,1,7,4,sample);           % \keppa_Th
S(:,6)=Beta_distribution(0.1,1,2,4,sample);           % \keppa_Tr
S(:,7)=Beta_distribution(0.1,1,3,4,sample);           % \keppa_Tc
S(:,8)=Beta_distribution(0.2,0.45,8,0.9,sample);       % \beta_Th
S(:,9)=Beta_distribution(0.1,0.6,0.1,0.5,sample);       % \beta_Tr
S(:,10)=Beta_distribution(0.1,1,5,7,sample);         % \beta_Tc

end