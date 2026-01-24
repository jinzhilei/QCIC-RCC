
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This function is used to describe the pharmacokinetic model of atezolizumab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A=PK_model_atezolizumab(~,x,V_1_A,V_2_A,Q_A,CI_A,V_max,K_A)

m = length(x); A = zeros(m,1);

A_1 = x(1); A_2 = x(2);

A(1) = ( 1 / V_1_A ) * ( Q_A * ( A_2 - A_1 ) - CI_A * A_1 - V_max * ( A_1 / ( A_1 + K_A ) ) );

A(2) = ( 1 / V_2_A ) * ( Q_A * ( A_1 - A_2 ));

end
