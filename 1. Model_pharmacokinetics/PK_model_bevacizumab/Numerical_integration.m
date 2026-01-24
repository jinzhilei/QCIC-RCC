
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This function is used to compute numerical integrals according to the trapezoidal rule

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = Numerical_integration(t,C)
    n = length(t);
    result = zeros(1,n-1);
    AUC = 0;
    for i = 1:n-1
        AUC = 0.5 * ( ( C(i) + C(i+1) ) * ( t(2) - t(1) ) )+AUC ; 
        result(i) = AUC;
    end
end
