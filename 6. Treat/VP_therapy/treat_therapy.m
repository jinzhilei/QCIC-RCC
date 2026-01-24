function treat_therapy(a,Drug_day,b,d,i,r,V_1_A,V_2_A,Q_A,CI_A,V_max_A,K_A_A,V_1_B,V_2_B,Q_B,CI_B,V_max_B,K_B_A,...
    theta_TN4,theta_TN8,theta_D0,theta_M,kappa_Th,kappa_Tr,kappa_Tc,beta_Th,beta_Tr,beta_Tc,G_X1,beta_X1,d_X1,G_X2,beta_X2,...
    d_X2,Lambda_A,rho_X2X1,K_A,G_X3,G_X4,beta_X3,beta_X4,d_X3,d_X4,Lambda_B,rho_X3X1,rho_X4X1,rho_X4X2,rho_X4X3,K_B)
t_start = a;
t_tau = Drug_day;
t_end = t_start+t_tau;
t_calibration = d;
m=ceil(t_calibration/t_tau);
Q_Ate_Bev = zeros(1,32);
Q_Ate_start = [0 0];
Q_Bev_start = [0 0];
x0=initialzation_VP_AB(Q_Ate_start,Q_Bev_start,r(i));
Drug = PSO(t_start,t_end,x0,Q_Ate_Bev,b,V_1_A,V_2_A,Q_A,CI_A,V_max_A,K_A_A,V_1_B,V_2_B,Q_B,CI_B,V_max_B,K_B_A,...
    theta_TN4,theta_TN8,theta_D0,theta_M,kappa_Th,kappa_Tr,kappa_Tc,beta_Th,beta_Tr,beta_Tc,G_X1,beta_X1,d_X1,G_X2,beta_X2,...
    d_X2,Lambda_A,rho_X2X1,K_A,G_X3,G_X4,beta_X3,beta_X4,d_X3,d_X4,Lambda_B,rho_X3X1,rho_X4X1,rho_X4X2,rho_X4X3,K_B);

Q_Ate = [Drug(1) 0];
Q_Bev = [Drug(2) 0];

Q_Ate_Bev(29) = Drug(1);
Q_Ate_Bev(31) = Drug(2);
x0=initialzation_VP_AB(Q_Ate,Q_Bev,r(i));

tout = zeros(1,d);
yout = zeros(d,32);
for j = 1:m
    if j < m
        drug(j,:) = Drug;
        [t,y] = ode45(@(t,x) QCIC_VP_AB(t,x,V_1_A,V_2_A,Q_A,CI_A,V_max_A,K_A_A,V_1_B,V_2_B,Q_B,CI_B,V_max_B,K_B_A,...
            theta_TN4,theta_TN8,theta_D0,theta_M,kappa_Th,kappa_Tr,kappa_Tc,beta_Th,beta_Tr,beta_Tc,G_X1,beta_X1,d_X1,G_X2,beta_X2,...
            d_X2,Lambda_A,rho_X2X1,K_A,G_X3,G_X4,beta_X3,beta_X4,d_X3,d_X4,Lambda_B,rho_X3X1,rho_X4X1,rho_X4X2,rho_X4X3,K_B),t_start:1:t_end, x0);
        nt = length(t);
        tout(t_start:1:t_end) = t_start:1:t_end;
        yout(t_start:1:t_end,:) = y(1:nt,:);
        t_start = t_end;
        t_end = t_start + t_tau;
        Drug = PSO(t_start,t_end,y(nt,:),Q_Ate_Bev,b,V_1_A,V_2_A,Q_A,CI_A,V_max_A,K_A_A,V_1_B,V_2_B,Q_B,CI_B,V_max_B,K_B_A,...
            theta_TN4,theta_TN8,theta_D0,theta_M,kappa_Th,kappa_Tr,kappa_Tc,beta_Th,beta_Tr,beta_Tc,G_X1,beta_X1,d_X1,G_X2,beta_X2,...
            d_X2,Lambda_A,rho_X2X1,K_A,G_X3,G_X4,beta_X3,beta_X4,d_X3,d_X4,Lambda_B,rho_X3X1,rho_X4X1,rho_X4X2,rho_X4X3,K_B);
        Q_Ate_Bev(29) = Drug(1);
        Q_Ate_Bev(31) = Drug(2);
        x0 = y(nt,:)+Q_Ate_Bev;
    else
        drug(j,:) = Drug;
        [t,y] = ode45(@(t,x) QCIC_VP_AB(t,x,V_1_A,V_2_A,Q_A,CI_A,V_max_A,K_A_A,V_1_B,V_2_B,Q_B,CI_B,V_max_B,K_B_A,...
            theta_TN4,theta_TN8,theta_D0,theta_M,kappa_Th,kappa_Tr,kappa_Tc,beta_Th,beta_Tr,beta_Tc,G_X1,beta_X1,d_X1,G_X2,beta_X2,...
            d_X2,Lambda_A,rho_X2X1,K_A,G_X3,G_X4,beta_X3,beta_X4,d_X3,d_X4,Lambda_B,rho_X3X1,rho_X4X1,rho_X4X2,rho_X4X3,K_B),t_start:t_calibration, x0);
        nt = length(t);
        tout(t_start:1:t_calibration) = t_start:1:t_calibration;
        yout(t_start:1:t_calibration,:) = y(1:nt,:);

        file=['./Adaptive_therapy/Adaptive_therapy_',num2str(Drug_day),'/',num2str(b),'/VPsample_AB',num2str(i),'.mat'];
        save(file,'yout');
        file=['./Drug/Drug_',num2str(Drug_day),'/',num2str(b),'/VPsample_AB',num2str(i),'.mat'];
        save(file,'drug');

    end
end
end