
function S = tumor_parameter(z,p,P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%This script file is used to perform parameter extraction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = zeros(z,8);
S(:,1) = Beta_distribution(5e11,5e13,p(1),p(2),z);   
S(:,2) = Beta_distribution(0.03,0.15,p(3),p(4),z);    
S(:,3) = Beta_distribution(0.05,0.2,p(5),p(6),z);       

S(:,4) = Beta_distribution(4e11,4e13,P(1),P(2),z);     
S(:,5) = Beta_distribution(0.01,0.18,P(3),P(4),z);     
S(:,6) = Beta_distribution(0.05,0.18,P(5),P(6),z);    
S(:,7) = Beta_distribution(0.001,0.5,P(7),P(8),z);        
S(:,8) = Beta_distribution(0.01,1,P(9),P(10),z);      

end

function y = Beta_distribution(a,b,c,d,e)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%This script file is used to perform beta distribution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = a + (b-a) * betarnd(c, d, [1, e]);
end
