function y = solve_one_patient(i,Par, S, Par_A,dose_A,t_start, t_calibration,t_end)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This script file is used to perform parallel computation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

check_interval = 21;
Q_Ate = zeros(1,30);
Q_Ate(29) = dose_A;
if t_calibration == t_end
    check_times = check_interval:check_interval:t_calibration;
else
    check_times = [check_interval:check_interval:t_end,t_end];
end

t_all = [];
y_all = [];
reset_times = [84, 126:42:t_end,t_end];

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

G_X2 =  S(i,4);
d_X2 =  S(i,5);
beta_X2 = S(i,6);
Lambda_A =  S(i,7);
rho_X2X1 =  S(i,8);

V_1_A       =  Par_A(i,1);
V_2_A       =  Par_A(i,2);
Q_A           =  Par_A(i,3);
CI_A          =  Par_A(i,4);
V_max_A  =  Par_A(i,5);
K_A_A      =  Par_A(i,6);
Dose = zeros(length(check_times),1);
Dose(1) = (dose_A);

x0 = initialzation_VP_A(i,dose_A);
for i =1:length(check_times)

    [t, y] = ode45(@(t,x) QCIC_VP_A(t, x, ...
        theta_TN4, theta_TN8, theta_D0, theta_M, ...
        beta_Th, beta_Tr, beta_Tc, ...
        G_X1, d_X1,beta_X1,...
        G_X2,d_X2,beta_X2,Lambda_A,rho_X2X1, ...
        V_1_A,V_2_A,Q_A,CI_A,V_max_A,K_A_A), ...
        [t_start,check_times(i)], x0);
    nt = length(t);
    t_all(t_start:1:check_times(i)) = t_start:1:check_times(i);
    y_all(t_start:1:check_times(i),:) = y(1:nt,:);
    t_start = t(end);

    if ismember(check_times(i), reset_times)
        if (((3 .* (y(end,25)+y(end,26)) .* 2.572e-9 .* 0.8) ./ (0.37 * 4 * 3.14)).^(1/3)-...
            ((3 .* y_all(1,25) .* 2.572e-9 .* 0.8) ./ (0.37 * 4 * 3.14)).^(1/3))/...
            ((3 .* y_all(1,25) .* 2.572e-9 .* 0.8) ./ (0.37 * 4 * 3.14)).^(1/3)>0.2
        Dose(i) = 0;
        Q_Ate(29) = Dose(i);
        else
            Dose(i) = dose_A;
            Q_Ate(29) = Dose(i);
        end
    else
        if i == 1
            Dose(i) = dose_A;
        else
        Dose(i) = Dose(i-1);
        end
        Q_Ate(29) = Dose(i);
    end

      x0 = Q_Ate+y(end,:);

end
y = y_all;


end
