function main_parallel(i,t_start,t_end,theta_TN4,theta_TN8,theta_D0,theta_M,kappa_Th,kappa_Tr,kappa_Tc,beta_Th,beta_Tr,beta_Tc)

x0 = initialzation_VP();

[~, x] = ode45(@(t,x) QCIC_VP(t,x,theta_TN4,theta_TN8,theta_D0,theta_M,kappa_Th,kappa_Tr,kappa_Tc,beta_Th,beta_Tr,beta_Tc),t_start:1:t_end,x0);

file=['./VP_parallel/','VP',num2str(i),'.mat'];

save(file,'x');

end