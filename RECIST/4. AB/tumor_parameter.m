function S = tumor_parameter(z,p_X1,P_X2,P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%This script file is used to perform parameter extraction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = zeros(z,8);
S(:,1) = Beta_distribution(5e11,5e13,p_X1(1),p_X1(2),z);       % G_X1
S(:,2) = Beta_distribution(0.03,0.15,p_X1(3),p_X1(4),z);      %  d_X1
S(:,3) = Beta_distribution(0.05,0.2,p_X1(5),p_X1(6),z);        % beta_X1

S(:,4) = Beta_distribution(4e11,4e13,P_X2(1),P_X2(2),z);       % G_X2      
S(:,5) = Beta_distribution(0.01,0.18,P_X2(3),P_X2(4),z);     % \d_X2  
S(:,6) = Beta_distribution(0.05,0.18,P_X2(5),P_X2(6),z);     % beta_X2  
S(:,7) = Beta_distribution(0.001,0.5,P_X2(7),P_X2(8),z);         % \Lambda_A   
S(:,8) = Beta_distribution(0.01,1,P_X2(9),P_X2(10),z);      % \rho_X2X1   


S(:,9) = Beta_distribution(4e11,4e13,P(1),P(2),z);      % G_X3         
S(:,10) = Beta_distribution(8e10,5e12,P(3),P(4),z);       % G_X4        
S(:,11) = Beta_distribution(0.05,0.18,P(5),P(6),z);      % \beta_X3    
S(:,12) = Beta_distribution(0.03,0.145,P(7),P(8),z);     % \beta_X4     
S(:,13) = Beta_distribution(0.018,0.16,P(9),P(10),z);     % \d_X3       
S(:,14) = Beta_distribution(0.02,0.15,P(11),P(12),z);     % \d_X4       
S(:,15) = Beta_distribution(0.01,50,P(13),P(14),z);     % \Lambda_B   
S(:,16) = Beta_distribution(0.01,1,P(15),P(16),z);      % \rho_X3X1   
S(:,17) = Beta_distribution(0.01,0.5,P(17),P(18),z);      % \rho_X4X1   
S(:,18) = Beta_distribution(0.01,0.8,P(19),P(20),z);      % \rho_X4X2   
S(:,19) = Beta_distribution(0.01,0.8,P(21),P(22),z);      % \rho_X4X3   


end


function y = Beta_distribution(a,b,c,d,e)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%This script file is used to perform beta distribution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = a + (b-a) * betarnd(c, d, [1, e]);
end