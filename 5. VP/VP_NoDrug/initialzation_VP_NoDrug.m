function y0=initialzation_VP_NoDrug(r) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Variable initial value (cell)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file='./VP_patient/VP_sample_immune.mat';
A=cell2mat(struct2cell(load(file)));

global par

par.D0_A_0 = A(1,1);          par.TN4_A_0 = A(1,7);     par.TN8_A_0 = A(1,16);     par.M_A_0 = A(1,22);                   

par.D0_B_0 = A(1,2);          par.TN4_B_0 = A(1,8);     par.Th_B_0 = A(1,10);       par.Tr_B_0 = A(1,13); 
par.TN8_B_0 = A(1,17);      par.Tc_B_0 = A(1,19);      par.M_B_0 = A(1,23);

par.D_C_0 = A(1,4);            par.TN4_C_0 = A(1,9);      par.Th_C_0 = A(1,11);          par.Tr_C_0 = A(1,14);
par.TN8_C_0 = A(1,18);     par.Tc_C_0 = A(1,20);  

par.D0_D_0 = A(1,3);         par.D_D_0 = A(1,5);           par.Th_D_0 = A(1,12);      par.Tr_D_0 = A(1,15);
par.Tc_D_0 = A(1,21);        par.TAM_D_0 = A(1,24);    par.X1_D_0 = r;           par.X2_D_0 = 0;
par.X3_D_0 = 0;              par.X4_D_0 = 0;              

par.D_E_0 = A(6); 

y0 = [par.D0_A_0;par.D0_B_0;par.D0_D_0;par.D_C_0;par.D_D_0;par.D_E_0;...
    par.TN4_A_0;par.TN4_B_0;par.TN4_C_0;par.Th_B_0;par.Th_C_0;par.Th_D_0;...
    par.Tr_B_0;par.Tr_C_0;par.Tr_D_0;par.TN8_A_0;par.TN8_B_0;par.TN8_C_0;...
    par.Tc_B_0;par.Tc_C_0;par.Tc_D_0;par.M_A_0;par.M_B_0;par.TAM_D_0;...
    par.X1_D_0;par.X2_D_0;par.X3_D_0;par.X4_D_0];
   
end