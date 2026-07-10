function y = solve_one_patient(i,Par, S, t_start, t_calibration)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This script file is used to perform parallel computation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    theta_TN4 = Par(i,1);
    theta_TN8 = Par(i,2);
    theta_D0 = Par(i,3);
    theta_M = Par(i,4);
    beta_Th = Par(i,5);
    beta_Tr = Par(i,6);
    beta_Tc = Par(i,7);
    
    G_X1 = S(i,1);
    d_X1 = S(i,2);
    beta_X1 = S(i,3);
    x0 = initialzation_VP_NoDrug(i);
    [~, y] = ode45(@(t,x) QCIC_VP_NoDrug(t, x, ...
        theta_TN4, theta_TN8, theta_D0, theta_M, ...
        beta_Th, beta_Tr, beta_Tc, ...
        G_X1, d_X1,beta_X1), ...
        t_start:1:t_calibration, x0);

end
