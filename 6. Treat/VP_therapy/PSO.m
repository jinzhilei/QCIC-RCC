function gbest = PSO(t_start,t_end,y0,Q_Ate_Bev,m,V_1_A,V_2_A,Q_A,CI_A,V_max_A,K_A_A,V_1_B,V_2_B,Q_B,CI_B,V_max_B,K_B_A,...
    theta_TN4,theta_TN8,theta_D0,theta_M,kappa_Th,kappa_Tr,kappa_Tc,beta_Th,beta_Tr,beta_Tc,G_X1,beta_X1,d_X1,G_X2,beta_X2,...
    d_X2,Lambda_A,rho_X2X1,K_A,G_X3,G_X4,beta_X3,beta_X4,d_X3,d_X4,Lambda_B,rho_X3X1,rho_X4X1,rho_X4X2,rho_X4X3,K_B)

n = 10; 
narvs = 2; 
c1 = 2;  
c2 = 2;  
w_max = 0.9;  
w_min = 0.4; 
K = 50; 
vmax = [0.4 0.4]; 
Drug_min = [0 0]; 
% Drug_max = [400 350;600 525;800 700]; % 
% Drug_max = [600 525]; % 
% Drug_max = [800 700]; % 
Drug_max = [1200 1050;1800 1575;2400 2100]; %
x = zeros(n,narvs);
for i = 1: narvs
    x(:,i) = round(Drug_min(i) + (Drug_max(m,i)-Drug_min(i))*rand(n,1));   
end
v = repmat(-vmax,n,1) + rand(n,narvs)*diag(2*vmax) ;  
fit = zeros(n,1); 
for i = 1:n  
    Q_Ate_Bev(29) = x(i,1);
    Q_Ate_Bev(31) = x(i,2);
    x0 = y0+Q_Ate_Bev;
    [t,y] = ode45(@(t,x) QCIC_VP_AB(t,x,V_1_A,V_2_A,Q_A,CI_A,V_max_A,K_A_A,V_1_B,V_2_B,Q_B,CI_B,V_max_B,K_B_A,...
        theta_TN4,theta_TN8,theta_D0,theta_M,kappa_Th,kappa_Tr,kappa_Tc,beta_Th,beta_Tr,beta_Tc,G_X1,beta_X1,d_X1,G_X2,beta_X2,...
        d_X2,Lambda_A,rho_X2X1,K_A,G_X3,G_X4,beta_X3,beta_X4,d_X3,d_X4,Lambda_B,rho_X3X1,rho_X4X1,rho_X4X2,rho_X4X3,K_B),t_start:1:t_end, x0);
    nt = length(t);
    fit(i) =  1*(y(nt,25)-y(1,25))+2*(y(nt,26)-y(1,26))+2*(y(nt,27)-y(1,27))+4*(y(nt,28)-y(1,28));
end
pbest = x;   
ffit = fit;
ind = find(fit == min(fit), 1);  
gbest = x(ind,:);  
git = fit(ind);

for d = 1:K  
    for i = 1:n   
        f_i = fit(i);  
        f_avg = sum(fit)/n;  
        f_min = min(fit);
        if f_i <= f_avg
            if f_avg ~= f_min  
                w = w_min + (w_max - w_min)*(f_i - f_min)/(f_avg - f_min);
            else
                w = w_min;
            end
        else
            w = w_max;
        end
        v(i,:) = round(w*v(i,:) + c1*rand(1)*(pbest(i,:) - x(i,:)) + c2*rand(1)*(gbest - x(i,:)));  
        for j = 1: narvs
            if v(i,j) < -vmax(j)
                v(i,j) = -vmax(j);
            elseif v(i,j) > vmax(j)
                v(i,j) = vmax(j);
            end
        end
        x(i,:) = x(i,:) + v(i,:); 
        for j = 1: narvs
            if x(i,j) < Drug_min(j)
                x(i,j) = Drug_min(j);
            elseif x(i,j) > Drug_max(m,j)
                x(i,j) = Drug_max(m,j);
            end
        end

        Q_Ate_Bev(29) = x(i,1);
        Q_Ate_Bev(31) = x(i,2);

        x0 = y0+Q_Ate_Bev;
        [~,y] = ode45(@(t,x) QCIC_VP_AB(t,x,V_1_A,V_2_A,Q_A,CI_A,V_max_A,K_A_A,V_1_B,V_2_B,Q_B,CI_B,V_max_B,K_B_A,...
            theta_TN4,theta_TN8,theta_D0,theta_M,kappa_Th,kappa_Tr,kappa_Tc,beta_Th,beta_Tr,beta_Tc,G_X1,beta_X1,d_X1,G_X2,beta_X2,...
            d_X2,Lambda_A,rho_X2X1,K_A,G_X3,G_X4,beta_X3,beta_X4,d_X3,d_X4,Lambda_B,rho_X3X1,rho_X4X1,rho_X4X2,rho_X4X3,K_B),t_start:1:t_end, x0);

        fit(i) =  1*(y(nt,25)-y(1,25))+2*(y(nt,26)-y(1,26))+2*(y(nt,27)-y(1,27))+4*(y(nt,28)-y(1,28));
        if fit(i) < ffit(i,:)   
            pbest(i,:) = x(i,:);  
        end
        if  fit(i) < git  
            gbest = pbest(i,:); 
        end
    end

end

end


