function y0=initialzation_basic() 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Variable initial value (cell)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global par

par.D0_A_0 = 5e7;      par.TN4_A_0 = 3.5e8;     par.TN8_A_0 = 3.5e8;     par.M_A_0 = 5e7;                   

par.D0_B_0 = 5e7;      par.TN4_B_0 = 3.5e8;     par.Th_B_0 = 5e7;        par.Tr_B_0 = 5e7; 
par.TN8_B_0 = 3.5e8;   par.Tc_B_0 = 5e7;        par.M_B_0 = 5e7;

par.D_C_0 = 5e7;       par.TN4_C_0 = 3.5e8;     par.Th_C_0 = 5e7;        par.Tr_C_0 = 5e7;
par.TN8_C_0 = 3.5e8;   par.Tc_C_0 = 5e7;  

par.D0_D_0 = 5e7;      par.D_D_0 = 5e7;         par.Th_D_0 = 5e7;        par.Tr_D_0 = 5e7;
par.Tc_D_0 = 5e7;      par.TAM_D_0 = 5e7;       par.X1_D_0 = 0;        par.X2_D_0 = 0;
par.X3_D_0 = 0;        par.X4_D_0 = 0;

par.D_E_0 = 5e7; 

y0 = [par.D0_A_0;par.D0_B_0;par.D0_D_0;par.D_C_0;par.D_D_0;par.D_E_0;...
    par.TN4_A_0;par.TN4_B_0;par.TN4_C_0;par.Th_B_0;par.Th_C_0;par.Th_D_0;...
    par.Tr_B_0;par.Tr_C_0;par.Tr_D_0;par.TN8_A_0;par.TN8_B_0;par.TN8_C_0;...
    par.Tc_B_0;par.Tc_C_0;par.Tc_D_0;par.M_A_0;par.M_B_0;par.TAM_D_0;
    par.X1_D_0;par.X2_D_0;par.X3_D_0;par.X4_D_0];
   
end