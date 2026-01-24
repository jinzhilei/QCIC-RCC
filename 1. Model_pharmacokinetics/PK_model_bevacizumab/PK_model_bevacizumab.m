
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This function is used to describe the pharmacokinetic model of bevacizumab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function B=PK_model_bevacizumab(~,x,V_1_B,V_2_B,Q_B,CI_B,V_max,K_B)

m = length(x); B = zeros(m,1);

B_1 = x(1); B_2 = x(2);

B(1) = ( 1 / V_1_B ) * ( Q_B * ( B_2 - B_1 ) - CI_B * B_1 - V_max * ( B_1 / ( B_1 + K_B ) ) );

B(2) = ( 1 / V_2_B ) * ( Q_B * ( B_1 - B_2 ));

end
