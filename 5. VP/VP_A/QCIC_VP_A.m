function y = QCIC_VP_A(t,x,theta_TN4,theta_TN8,theta_D0,theta_M,kappa_Th,kappa_Tr,kappa_Tc,beta_Th,beta_Tr,beta_Tc,G_X1,beta_X1,d_X1,G_X2,beta_X2,d_X2,...
                         Lambda_A,rho_X2X1,K_A, V_1_A,V_2_A,Q_A,CI_A,V_max_A,K_A_A)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% QCIC model (basic)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set: parameter
global par 
parameter_VP_A();

%% Prestore: variable
m=length(x);
y=zeros(m,1);

%% Variable substitution
D0_A = x(1);    D0_B = x(2);    D0_D = x(3);    D_C = x(4);    D_D = x(5);      D_E = x(6);
TN4_A = x(7);   TN4_B = x(8);   TN4_C = x(9);   Th_B = x(10);  Th_C = x(11);    Th_D = x(12);
Tr_B = x(13);   Tr_C = x(14);   Tr_D = x(15);   TN8_A = x(16); TN8_B = x(17);   TN8_C = x(18);
Tc_B = x(19);   Tc_C = x(20);   Tc_D = x(21);   M_A = x(22);   M_B = x(23);     TAM_D = x(24);
X1_D = x(25);   X2_D = x(26);   X3_D = x(27);   X4_D = x(28); A_B = x(29);     A_D = x(30);

%% Load Protein module (fast time scale) 
X_D = X1_D + X2_D + X3_D +X4_D;
% Cytokines and chemokines in compartment C

I2_C = (1 / par.d_I2) * ( par.delta_I2Th * Th_C); % IL2 in compartment C
I10_C = (1 / par.d_I10) * ( par.delta_I10Tr * Tr_C); % IL10 in compartment C
I12_C = (1 / par.d_I12) * ( par.delta_I12D * D_C); % IL12 in compartment C
Igamma_C = (1 / par.d_Igamma) * ( par.delta_IgammaTh * Th_C + par.delta_IgammaTc * Tc_C); % Igamma in compartment C
Tbeta_C = (1 / par.d_Tbeta) * ( par.delta_TbetaTr * Tr_C ); % Tbeta in compartment C
L19_C = (1 / par.d_L19) * ( par.delta_L19D * D_C ); % CCL19 in compartment C
L21_C = (1 / par.d_L21) * ( par.delta_L21D * D_C ); % CCL21 in compartment C

% Chemokines in compartment D

C9_D = (1 / par.d_C9) * ( par.delta_C9D * D_D + par.delta_C9T * Tc_D + par.delta_C9T * Th_D ); % CXCL9 in compartment D
C10_D = (1 / par.d_C10) * ( par.delta_C10D * D_D + par.delta_C10T * Tc_D + par.delta_C10T * Th_D ); % CXCL10 in compartment D
L2_D = (1 / par.d_L2) * ( par.delta_L2X * X_D ); % CCL2 in compartment D
L17_D = (1 / par.d_L17) * ( par.delta_L17TAM * TAM_D ); % CCL17 in compartment D
L20_D = (1 / par.d_L20) * ( par.delta_L20X * X_D ); % CCL20 in compartment D
L22_D = (1 / par.d_L22) * ( par.delta_L22TAM * TAM_D); % CCL22 in compartment D

%% Load Protein module (slow time scale)

% ====== DSTC ======

proliferation_X1_D = (beta_X1 / (1 + par.Lambda_B * (par.B_D^par.n_B / (par.P^par.n_B+par.B_D^par.n_B) ) ) )* (1 - (par.a11 * X1_D + par.a12 * X2_D + par.a13 * X3_D + par.a14 * X4_D)/G_X1) * X1_D;
transformation_X2X1 = rho_X2X1 * (A_D^par.n_3 / (K_A^par.n_3 + A_D^par.n_3) ) * X1_D;
transformation_X3X1 = par.rho_X3X1 * ( par.B_D^par.n_3 / (par.K_B^par.n_3 + par.B_D^par.n_3) ) * X1_D;
transformation_X4X1 = par.rho_X4X1 * ( A_D^par.n_3 / (K_A^par.n_3 + A_D^par.n_3) )*( par.B_D^par.n_3 / (par.K_B^par.n_3 + par.B_D^par.n_3) ) * X1_D;
apoptosis_X1_D = d_X1 * X1_D;
immune_X1_D = (par. eta_Th * Th_D + par.eta_Tc * Tc_D) * ( 1 + Lambda_A * (A_D^par.n_A / (par.Q^par.n_A+A_D^par.n_A) ) ) * (1/ (Tr_D/par.K_Tr + 1)) * X1_D;

% ====== ARTC ======
proliferation_X2_D = (beta_X2 / (1 + par.Lambda_B * (par.B_D^par.n_B / (par.P^par.n_B+par.B_D^par.n_B) ) ) ) * (1 - (par.a21 * X1_D + par.a22 * X2_D + par.a23 * X3_D +par.a24 * X4_D)/G_X2) * X2_D;
transformation_X4X2 = par.rho_X4X2 * ( par.B_D^par.n_4 / (par.K_B^par.n_4 + par.B_D^par.n_4) ) * X2_D;
apoptosis_X2_D = d_X2 * X2_D;
immune_X2_D = (par.eta_Th * Th_D + par.eta_Tc * Tc_D) * (  (1 / (Tr_D/par.K_Tr + 1))) * X2_D;

% ====== BRTC ======
proliferation_X3_D = par.beta_X3 * (1 - (par.a31 * X1_D + par.a32 * X2_D + par.a33 * X3_D + par.a34 * X4_D )/par.G_X3) * X3_D;
transformation_X4X3 = par.rho_X4X3 * ( A_D^par.n_3 / (K_A^par.n_3 + A_D^par.n_3) ) * X3_D;
apoptosis_X3_D = par.d_X3 * X3_D;
immune_X3_D = (par. eta_Th * Th_D + par.eta_Tc * Tc_D) *  ( 1 + Lambda_A * (A_D^par.n_A / (par.Q^par.n_A+A_D^par.n_A) ) ) * (1/ (Tr_D/par.K_Tr + 1)) * X3_D;

% ====== DRTC ======
proliferation_X4_D = par.beta_X4 * (1 - (par.a41 * X1_D + par.a42 * X2_D + par.a43 * X3_D + par.a44 * X4_D )/par.G_X4) * X4_D;
apoptosis_X4_D = par.d_X4 * X4_D;
immune_X4_D = (par.eta_Th * Th_D + par.eta_Tc * Tc_D) * (  (1 / (Tr_D/par.K_Tr + 1))) * X4_D;

% ====== TAA ======
TA_D = (1 / par.d_TA) * ((par.u_d * apoptosis_X1_D + par.u_d * apoptosis_X2_D + par.u_d * apoptosis_X3_D + par.u_d * apoptosis_X4_D) + (par.u_I * immune_X1_D + par.u_I * immune_X2_D + par.u_I * immune_X3_D + par.u_I * immune_X4_D));

% ====== D_0 ======
source_D0_A = theta_D0;
migration_D0_A = par.V_D0M_BA * (par.V_B / par.V_A) * D0_B - par.V_D0M_AB * D0_A;
apoptosis_D0_A = par.d_D0 * D0_A;

migration_D0_B = par.V_D0M_AB * (par.V_A / par.V_B) * D0_A - par.V_D0M_BA* D0_B + par.V_D0M_DB * (par.V_D / par.V_B) * D0_D - par.V_D0M_BD * D0_B;
chemotaxis_D0_BD = par.X_D0L20_BD * (L20_D / (par.K_L20 + L20_D)) * D0_B;
apoptosis_D0_B = par.d_D0 * D0_B;

migration_D0_D = par.V_D0M_BD * (par.V_B / par.V_D) * D0_B - par.V_D0M_DB * D0_D;
tranformation_D0_D = par.rho_DD0 * (TA_D^par.n_1 / (par.K_TA^par.n_1 + TA_D^par.n_1 )) * D0_D;
apoptosis_D0_D = par.d_D0 * D0_D;

% ====== D ======

migration_D_D = par.V_D_DE * D_D;
apoptosis_D_D = par.d_D * D_D;

migration_D_E = par.V_D_DE * (par.V_D / par.V_E) * D_D - par.V_D_EC * D_E;
chemotaxis_D_EC = (par.X_DL19_EC * (L19_C / (par.K_L19 + L19_C)) + par.X_DL21_EC * (L21_C / (par.K_L21 + L21_C))) * D_E;
apoptosis_D_E = par.d_D * D_E;

migration_D_C = par.V_D_EC * (par.V_E / par.V_C) * D_E;
apoptosis_D_C = par.d_D * D_C;

% ====== TN4 ======

source_TN4_A = theta_TN4;
migration_TN4_A = par.V_TN_BA * (par.V_B / par.V_A) * TN4_B - par.V_TN_AB * TN4_A;
apoptosis_TN4_A = par.d_TN4 * TN4_A;

migration_TN4_B = (par.V_TN_AB * (par.V_A / par.V_B) * TN4_A - par.V_TN_BA * TN4_B) + (par.V_TN_CB * (par.V_C / par.V_B) * TN4_C - par.V_TN_BC * TN4_B);
chemotaxis_TN4_BC = (par.X_TL19_BC * (L19_C / (par.K_L19 + L19_C)) + par.X_TL21_BC * (L21_C / (par.K_L21 + L21_C))) * TN4_B;
apoptosis_TN4_B = par.d_TN4 * TN4_B;

proliferation_TN4_C = par.beta_TN4 * TN4_C * (1 - (TN4_C / par.G_TN4));
migration_TN4_C = par.V_TN_BC * (par.V_B / par.V_C) * TN4_B - par.V_TN_CB * TN4_C;
differentiation_Th_C = kappa_Th * (1 + par.lambda_ThI12 * (I12_C / (par.K_I12 + I12_C)) + par.lambda_ThIgamma * (Igamma_C / (par.K_Igamma + Igamma_C))) * (1 / (1 + I10_C / par.K_I10)) * (1 / (1 + Tr_C / par.K_Tr)) * (D_C^par.n_2 / (par.K_D^par.n_2 + D_C^par.n_2))  * TN4_C;
differentiation_Tr_C = kappa_Tr * (1 + par.lambda_TrTbeta * (Tbeta_C / (par.K_Tbeta + Tbeta_C))) * (1 / (1 + Igamma_C / par.K_TrIgamma)) * (D_C^par.n_2 / (par.K_D^par.n_2 + D_C^par.n_2))  * TN4_C ;
apoptosis_TN4_C = par.d_TN4 * TN4_C;

% ====== Th ======

migration_Th_B = (par.V_T_CB * (par.V_C / par.V_B) * Th_C - par.V_T_BC * Th_B) + (par.V_T_DB * (par.V_D / par.V_B) * Th_D - par.V_T_BD * Th_B);
chemotaxis_Th_BD = (par.X_ThC9_BD * (C9_D / (par.K_C9 + C9_D)) + par.X_ThC10_BD * (C10_D / (par.K_C10 + C10_D))) * Th_B;
apoptosis_Th_B = par.d_Th * Th_B;

proliferation_Th_C = beta_Th * (1 + par.lambda_I2 * (I2_C / (par.K_I2 + I2_C))) * Th_C ;
migration_Th_C = (par.V_T_BC * (par.V_B / par.V_C) * Th_B - par.V_T_CB * Th_C);
apoptosis_Th_C= par.d_Th * Th_C;

migration_Th_D = par.V_T_BD * (par.V_B / par.V_D) * Th_B - par.V_T_DB * Th_D;
apoptosis_Th_D = par.d_Th * Th_D;

% ====== Tr ======

migration_Tr_B = (par.V_T_CB * (par.V_C / par.V_B) * Tr_C - par.V_T_BC * Tr_B) + (par.V_T_DB * (par.V_D / par.V_B) * Tr_D - par.V_T_BD * Tr_B);
chemotaxis_Tr_BD = (par.X_TrL17_BD * (L17_D / (par.K_L17 + L17_D)) + par.X_TrL22_BD * (L22_D / (par.K_L22 + L22_D))) * Tr_B;
apoptosis_Tr_B = par.d_Tr * Tr_B;

proliferation_Tr_C = beta_Tr * Tr_C ;
migration_Tr_C = (par.V_T_BC * (par.V_B / par.V_C) * Tr_B - par.V_T_CB * Tr_C);
apoptosis_Tr_C= par.d_Tr * Tr_C;

migration_Tr_D = par.V_T_BD * (par.V_B / par.V_D) * Tr_B - par.V_T_DB * Tr_D;
apoptosis_Tr_D = par.d_Tr * Tr_D;

% ====== TN8 ======

source_TN8_A = theta_TN8;
migration_TN8_A = par.V_TN_BA * (par.V_B / par.V_A) * TN8_B - par.V_TN_AB * TN8_A;
apoptosis_TN8_A = par.d_TN8 * TN8_A;

migration_TN8_B = (par.V_TN_AB * (par.V_A / par.V_B) * TN8_A - par.V_TN_BA * TN8_B) + (par.V_TN_CB * (par.V_C / par.V_B) * TN8_C - par.V_TN_BC * TN8_B);
chemotaxis_TN8_BC = (par.X_TL19_BC * (L19_C / (par.K_L19 + L19_C)) + par.X_TL21_BC * (L21_C / (par.K_L21 + L21_C))) * TN8_B;
apoptosis_TN8_B = par.d_TN8 * TN8_B;

proliferation_TN8_C = par.beta_TN8 * TN8_C * (1 - (TN8_C / par.G_TN8));
migration_TN8_C = par.V_TN_BC * (par.V_B / par.V_C) * TN8_B - par.V_TN_CB * TN8_C;
differentiation_Tc_C = kappa_Tc * (1 + par.lambda_TcI12 * (I12_C / (par.K_I12 + I12_C))) * (1 / (1 + I10_C / par.K_I10)) * (1 / (1 + Tr_C / par.K_Tr)) * (D_C^par.n_2 / (par.K_D^par.n_2 + D_C^par.n_2)) * TN8_C;
apoptosis_TN8_C = par.d_TN8 * TN8_C;

% ====== Tc ======

migration_Tc_B = (par.V_T_CB * (par.V_C / par.V_B) * Tc_C - par.V_T_BC * Tc_B) + (par.V_T_DB * (par.V_D / par.V_B) * Tc_D - par.V_T_BD * Tc_B);
chemotaxis_Tc_BD = (par.X_TcC9_BD * (C9_D / (par.K_C9 + C9_D)) + par.X_TcC10_BD * (C10_D / (par.K_C10 + C10_D))) * Tc_B;
apoptosis_Tc_B = par.d_Tc * Tc_B;

proliferation_Tc_C = beta_Tc * (1 + par.lambda_I2 * (I2_C / (par.K_I2 + I2_C))) * Tr_C ;
migration_Tc_C = (par.V_T_BC * (par.V_B / par.V_C) * Tc_B - par.V_T_CB * Tc_C);
apoptosis_Tc_C= par.d_Tc * Tc_C;

migration_Tc_D = par.V_T_BD * (par.V_B / par.V_D) * Tc_B - par.V_T_DB * Tc_D;
apoptosis_Tc_D = par.d_Tc * Tc_D;

% ====== M ======
source_M_A = theta_M;
migration_M_A = (par.V_D0M_BA * (par.V_B / par.V_A) * M_B - par.V_D0M_AB * M_A);
apoptosis_M_A = par.d_M * M_A;

migration_M_B = (par.V_D0M_AB * (par.V_A / par.V_B) * M_A - par.V_D0M_BA * M_B) + (par.V_D0M_DB * (par.V_D / par.V_B) * TAM_D - par.V_D0M_BD * M_B);
chemotaxis_M_BD = (par.X_ML2_BD * (L2_D / (par.K_L2 + L2_D))) * M_B;
apoptosis_M_B = par.d_M * M_B;

migration_TAM_D = (par.V_D0M_BD * (par.V_B / par.V_D) * M_B - par.V_D0M_DB * TAM_D);
apoptosis_TAM_D = par.d_TAM * TAM_D;

%% ODE equations combination

% ODE of D0

y(1) = source_D0_A + migration_D0_A - apoptosis_D0_A;
y(2) = migration_D0_B - chemotaxis_D0_BD - apoptosis_D0_B;
y(3) = migration_D0_D + chemotaxis_D0_BD - tranformation_D0_D -apoptosis_D0_D;

% ODE of D

y(4) = migration_D_C + chemotaxis_D_EC - apoptosis_D_C;
y(5) = tranformation_D0_D - migration_D_D - apoptosis_D_D;
y(6) = migration_D_E - chemotaxis_D_EC - apoptosis_D_E;

% ODE of TN4

y(7) = source_TN4_A + migration_TN4_A - apoptosis_TN4_A;
y(8) = migration_TN4_B - chemotaxis_TN4_BC - apoptosis_TN4_B;
y(9) = proliferation_TN4_C + migration_TN4_C + chemotaxis_TN4_BC - differentiation_Th_C - differentiation_Tr_C - apoptosis_TN4_C;

% ODE of Th

y(10) = migration_Th_B - chemotaxis_Th_BD - apoptosis_Th_B;
y(11) = proliferation_Th_C + migration_Th_C + differentiation_Th_C - apoptosis_Th_C;
y(12) = migration_Th_D + chemotaxis_Th_BD - apoptosis_Th_D;

% ODE of Tr

y(13) = migration_Tr_B - chemotaxis_Tr_BD - apoptosis_Tr_B;
y(14) = proliferation_Tr_C + migration_Tr_C + differentiation_Tr_C -apoptosis_Tr_C;
y(15) = migration_Tr_D + chemotaxis_Tr_BD - apoptosis_Tr_D;

% ODE of TN8

y(16) = source_TN8_A + migration_TN8_A - apoptosis_TN8_A;
y(17) = migration_TN8_B - chemotaxis_TN8_BC - apoptosis_TN8_B;
y(18) = proliferation_TN8_C + migration_TN8_C + chemotaxis_TN8_BC - differentiation_Tc_C -apoptosis_TN8_C;

% ODE of Tc

y(19) = migration_Tc_B - chemotaxis_Tc_BD - apoptosis_Tc_B;
y(20) = proliferation_Tc_C + migration_Tc_C + differentiation_Tc_C - apoptosis_Tc_C;
y(21) = migration_Tc_D + chemotaxis_Tc_BD - apoptosis_Tc_D;

% ODE of M

y(22) = source_M_A + migration_M_A - apoptosis_M_A;
y(23) = migration_M_B - chemotaxis_M_BD - apoptosis_M_B;
y(24) = migration_TAM_D + chemotaxis_M_BD - apoptosis_TAM_D;

% ODE of tumor

y(25) = proliferation_X1_D - apoptosis_X1_D - transformation_X2X1 - transformation_X3X1 - transformation_X4X1 - immune_X1_D;
y(26) = proliferation_X2_D - apoptosis_X2_D + transformation_X2X1 - transformation_X4X2 - immune_X2_D;
y(27) = proliferation_X3_D - apoptosis_X3_D + transformation_X3X1 - transformation_X4X3 - immune_X3_D;
y(28) = proliferation_X4_D - apoptosis_X4_D + transformation_X4X1 + transformation_X4X2 + transformation_X4X3 - immune_X4_D;

%% atezolizumab

y(29) = ( 1 / V_1_A ) * ( Q_A * ( A_D - A_B ) - CI_A * A_B - V_max_A * ( A_B / ( A_B + K_A_A ) ) );
y(30) = ( 1 / V_2_A ) * ( Q_A * ( A_B - A_D ));

end