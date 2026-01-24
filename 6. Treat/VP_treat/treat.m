function treat(a,c,d,i,r,p,m,n,Drug_Ate_init,Drug_Bev_init,ratio,Drug,V_1_A,V_2_A,Q_A,CI_A,V_max_A,K_A_A,V_1_B,V_2_B,Q_B,CI_B,V_max_B,K_B_A,...
                        theta_TN4,theta_TN8,theta_D0,theta_M,kappa_Th,kappa_Tr,kappa_Tc,beta_Th,beta_Tr,beta_Tc,G_X1,beta_X1,d_X1,G_X2,beta_X2,...
                    d_X2,Lambda_A,rho_X2X1,K_A,G_X3,G_X4,beta_X3,beta_X4,d_X3,d_X4,Lambda_B,rho_X3X1,rho_X4X1,rho_X4X2,rho_X4X3,K_B)

                    t_start = a;
                    t_tau = c;
                    t_end = t_start+t_tau;
                    t_calibration = d;
                    o=ceil(t_calibration/t_tau);

                    Q_Ate_start = [Drug_Ate_init 0];
                    Q_Bev_start = [Drug_Bev_init 0];
                    Q_Ate_Bev = zeros(1,32);

                    x0=initialzation_VP_AB(Q_Ate_start,Q_Bev_start,r(i));
                    tout = zeros(1,d);
                    yout = zeros(d,32);
                    for j = 1:o
                        if j < o
                            [t,y] = ode45(@(t,x) QCIC_VP_AB(t,x,V_1_A,V_2_A,Q_A,CI_A,V_max_A,K_A_A,V_1_B,V_2_B,Q_B,CI_B,V_max_B,K_B_A,...
                        theta_TN4,theta_TN8,theta_D0,theta_M,kappa_Th,kappa_Tr,kappa_Tc,beta_Th,beta_Tr,beta_Tc,G_X1,beta_X1,d_X1,G_X2,beta_X2,...
                    d_X2,Lambda_A,rho_X2X1,K_A,G_X3,G_X4,beta_X3,beta_X4,d_X3,d_X4,Lambda_B,rho_X3X1,rho_X4X1,rho_X4X2,rho_X4X3,K_B),t_start:1:t_end, x0);
                            nt = length(t);
                            tout(t_start:1:t_end) = t_start:1:t_end;
                            yout(t_start:1:t_end,:) = y(1:nt,:);
                            t_start = t_end;
                            t_end = t_start + t_tau;
                            if mod(j+1,3) == 1
                                Q_Ate_Bev(29) = ratio(p,1) * Drug(1,m);
                                Q_Ate_Bev(31) = ratio(p,1) * Drug(2,m);
                            elseif mod(j+1,3) == 2
                                Q_Ate_Bev(29) = ratio(p,2) * Drug(1,m);
                                Q_Ate_Bev(31) = ratio(p,2) * Drug(2,m);
                            elseif mod(j+1,3) == 0
                                Q_Ate_Bev(29) = ratio(p,3) * Drug(1,m);
                                Q_Ate_Bev(31) = ratio(p,3) * Drug(2,m);
                            end
                            x0 = Q_Ate_Bev+y(end,:);
                        else
                            [t,y] = ode45(@(t,x) QCIC_VP_AB(t,x,V_1_A,V_2_A,Q_A,CI_A,V_max_A,K_A_A,V_1_B,V_2_B,Q_B,CI_B,V_max_B,K_B_A,...
                        theta_TN4,theta_TN8,theta_D0,theta_M,kappa_Th,kappa_Tr,kappa_Tc,beta_Th,beta_Tr,beta_Tc,G_X1,beta_X1,d_X1,G_X2,beta_X2,...
                    d_X2,Lambda_A,rho_X2X1,K_A,G_X3,G_X4,beta_X3,beta_X4,d_X3,d_X4,Lambda_B,rho_X3X1,rho_X4X1,rho_X4X2,rho_X4X3,K_B),t_start:t_calibration, x0);
                            nt = length(t);
                            tout(t_start:1:t_calibration) = t_start:1:t_calibration;
                            yout(t_start:1:t_calibration,:) = y(1:nt,:);
                        end
                    end
                    file=strcat('treat_AB/treat',num2str(p),'/(',num2str(m),',',num2str(n),')','/VPsample_AB',num2str(i),'.mat');
                    save(file,'yout');

end