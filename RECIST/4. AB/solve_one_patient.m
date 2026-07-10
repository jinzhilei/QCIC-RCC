function y = solve_one_patient(i,Par, S,Dose_par,dose_A,dose_B,t_start, t_calibration,t_end)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This script file is used to perform parallel computation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

check_interval = 21;
Q_dose = zeros(1,32);
Q_dose(29) = dose_A;
Q_dose(31) = dose_B;
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

G_X3 =  S(i,9);
G_X4 = S(i,10);
beta_X3 = S(i,11);
beta_X4 = S(i,12);
d_X3 =  S(i,13);
d_X4 =  S(i,14);
Lambda_B =  S(i,15);
rho_X3X1 =  S(i,16);
rho_X4X1 =  S(i,17);
rho_X4X2 =  S(i,18);
rho_X4X3 =  S(i,19);

V_1_A       =  Dose_par(i,1);
V_2_A       =  Dose_par(i,2);
Q_A           =  Dose_par(i,3);
CI_A          =  Dose_par(i,4);
V_max_A  =  Dose_par(i,5);
K_A_A      =  Dose_par(i,6);

V_1_B       =  Dose_par(i,7);
V_2_B       =  Dose_par(i,8);
Q_B           =  Dose_par(i,9);
CI_B          =  Dose_par(i,10);
V_max_B  =  Dose_par(i,11);
K_B_A      =  Dose_par(i,12);

Dose = zeros(length(check_times),2);
Dose(1,1) = (dose_A);
Dose(1,2) = (dose_B);

x0 = initialzation_VP_AB(i,dose_A,dose_B);
for i = 1:length(check_times)

    [t, y] = ode45(@(t,x) QCIC_VP_AB(t, x, ...
        theta_TN4, theta_TN8, theta_D0, theta_M, ...
        beta_Th, beta_Tr, beta_Tc, ...
        G_X1, d_X1,beta_X1,...
        G_X2,d_X2,beta_X2,Lambda_A,rho_X2X1, ...
        G_X3,G_X4,d_X3,d_X4,beta_X3,beta_X4,Lambda_B,rho_X3X1,rho_X4X1,rho_X4X2,rho_X4X3,...
        V_1_A,V_2_A,Q_A,CI_A,V_max_A,K_A_A,...
        V_1_B,V_2_B,Q_B,CI_B,V_max_B,K_B_A), ...
        t_start:check_times(i), x0);
    
    nt = length(t);
    t_all(t_start:1:check_times(i)) = t_start:1:check_times(i);
    y_all(t_start:1:check_times(i),:) = y(1:nt,:);
    t_start = t(end);

    if ismember(check_times(i), reset_times)
        if (((3 .* (y(end,25)+y(end,26)+y(end,27)+y(end,28)) .* 2.572e-9 .* 0.8) ./ (0.37 * 4 * 3.14)).^(1/3)-...
            ((3 .* y_all(1,25) .* 2.572e-9 .* 0.8) ./ (0.37 * 4 * 3.14)).^(1/3))/...
            ((3 .* y_all(1,25) .* 2.572e-9 .* 0.8) ./ (0.37 * 4 * 3.14)).^(1/3)<0.2 && any(Dose(i-1,:)~=0) 
            Dose(i,1) = dose_A;
            Dose(i,2) = dose_B;
            Q_dose(29) = Dose(i,1);
            Q_dose(31) = Dose(i,2);
            
        else
            Dose(i,1) = 0;
            Dose(i,2) = 0;
            Q_dose(29) = Dose(i,1);
            Q_dose(31) = Dose(i,2);
        end
    else
        if i == 1
            Dose(i,1) = dose_A;
            Dose(i,2) = dose_B;
        else
        Dose(i,:) = Dose(i-1,:);
        end
        Q_dose(29) = Dose(i,1);
        Q_dose(31) = Dose(i,2);
    end

      x0 = Q_dose+y(end,:);

end
y = y_all;


end
