function parameter_VP_NoDrug()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Parameter (VP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global par

%% kinetic parameters of migration 
par.V_D0M_BA = 0.05;
par.V_D0M_AB = 0.5;
par.V_D0M_BD = 0.1;
par.V_D0M_DB = 0.2;

par.V_D_DE = 0.01;
par.V_D_EC = 0.05;

par.V_TN_BA = 0.04;
par.V_TN_AB = 0.5;
par.V_TN_CB = 0.2;
par.V_TN_BC = 1e-3;

par.V_T_CB = 0.625;
par.V_T_BC = 1e-3;
par.V_T_DB = 2.00;
par.V_T_BD = 0.65;
%% kinetic parameters of volume   

par.V_A = 3;
par.V_B = 4.50;%
par.V_C = 0.01;%
par.V_D = 0.8;%
par.V_E = 0.1;

%% kinetic parameters of chemotaxis  

par.X_D0L20_BD = 0.02;           % 0.02 Sozzani.JImmunol.1998
par.X_DL19_EC = 0.1;
par.X_DL21_EC = 0.1;
par.X_TL19_BC = 0.1;
par.X_TL21_BC = 0.5;             % 0.5  Gunn.PNAS.1998

par.X_ThC9_BD = 0.09;            % 0.09 Loetscher.JEM.1996
par.X_ThC10_BD = 0.15;           % 0.15 Loetscher.JEM.1996
par.X_TrL17_BD  = 0.15;          % 0.15 Iellem.JEM.2001
par.X_TrL22_BD = 0.17;           % 0.17 Iellem.JEM.2001
par.X_TcC9_BD = 0.09;            % 0.09 Loetscher.JEM.1996
par.X_TcC10_BD = 0.12;           % 0.12 Loetscher.JEM.1996
par.X_ML2_BD = 0.35;             % 0.35 Bleul.JEM.1996

%% kinetic parameters of apoptosis  

par.d_D0 = 0.1;                 % 0.1 Lai.PNAS.2018
par.d_D = 0.1;                  % 0.1 Lai.PNAS.2018

par.d_M = 0.3;                  % 0.3 Li.Heliyon.2022
par.d_TAM = 0.1;                

par.d_TA = 0.01;                % 0.01 Zandarashvili.JBiolChem.2013
par.d_I2 = 2.38;                % 2.38 Lai.PNAS.2018
par.d_I10 = 8.32;               % 8.32 Lai.PNAS.2018
par.d_I12 = 1.38;               % 1.38 Lai.PNAS.2018
par.d_Tbeta = 499.07;           % 499.07 Lai.PNAS.2018
par.d_Igamma = 7.68;            % 7.68 Wang.JITC.2021
par.d_L2 = 1.73;                % 1.73 Wang.JITC.2021
par.d_L17 = 20.00;
par.d_L19 = 20.00;
par.d_L20 = 20.00;
par.d_L21 = 20.00;
par.d_L22 = 288.81;             % 288.81 Lambeir.JBiolChem.2001
par.d_C9 = 40.77;               % 40.77 Lambeir.JBiolChem.2001
par.d_C10 = 256.72;             % 256.72 Lambeir.JBiolChem.2001
%% kinetic parameters of half-saturation  

par.K_TA = 1e10;
par.K_I2 = 6e3;                 % 5e3 Li.FrontImmunol.2023
par.K_I12 = 8e3;                % 8e3 Lai.PNAS.2018
par.K_Igamma = 1e3;             % 1e3 Li.FrontImmunol.2023
par.K_Tbeta = 2.06e4;           % 2.06e4 Lai.SciChinaMath.2020

par.K_L2 = 7.5e4;               % 7.5e4 Li.FrontImmunol.2023
par.K_L17 = 2e4;              % 2.5e4 Li.FrontImmunol.2023
par.K_L19 = 2e5;                % 2e5 Li.FrontImmunol.2023
par.K_L20 = 8e3;                % 8e3 Li.FrontImmunol.2023
par.K_L21 = 2e7;                % 2e7 Li.FrontImmunol.2023
par.K_L22 = 5e5;                % 5e5 Li.FrontImmunol.2023
par.K_C9 = 5e4;                 % 5e4 Li.FrontImmunol.2023
par.K_C10 = 1.75e5;             % 1.75e5 Li.FrontImmunol.2023

par.K_D = 7e8;

par.K_I10 = 4.38e4;             % 4.38e4 Lai.PNAS.2018
par.K_TrIgamma = 1e4;
par.K_Tr = 5e8;

%% Regulation coefficient in the form of Hill or Michaelis-Menten function

par.lambda_ThIgamma = 0.6;
par.lambda_ThI12 = 0.6;
par.lambda_TrTbeta = 0.5;
par.lambda_TcI12 = 0.8;
par.lambda_I2 = 0.25;%

%% kinetic parameters of secretion (cytokine)  

par.delta_I2Th = 1e-4;           % 2e-4 Li.Heliyon.2022
par.delta_I10Tr = 2e-4;          % 1.4e-5 Robertson.JTB.2012
par.delta_I12D = 1e-4;
par.delta_TbetaTr = 9e-3;
par.delta_IgammaTh = 1e-6;
par.delta_IgammaTc = 1e-6;

par.delta_L2X = 1e-5;            % 1e-5 Lai.PNAS.2018
par.delta_L17TAM = 4e-3;         % 4e-4 Lai.PNAS.2018
par.delta_L19D = 4e-3;
par.delta_L20X = 3e-4;
par.delta_L21D = 0.2;
par.delta_L22TAM = 1;
par.delta_C9D = 2.5e-3;
par.delta_C9T = 2e-2;
par.delta_C10D = 2e-1;
par.delta_C10T = 1e-2;

%% kinetic parameters of antigen presentation  

par.n_1 = 2;
par.n_2 = 3;
par.n_3 = 4;
par.n_4 = 4;
par.n_A = 8;
par.n_B = 8;

par.a11 = 1;%
par.a12 = 0.75;%
par.a13 = 0.75;%
par.a14 = 0.75;
par.a21 = 1.25;%
par.a22 = 1;%
par.a23 = 1;%
par.a24 = 0.75;
par.a31 = 1.25;%
par.a32 = 1;%
par.a33 = 1;%
par.a34 = 0.75;
par.a41 = 1.25;
par.a42 = 1.25;
par.a43 = 1.25;
par.a44 = 1;

par.Q = 300;
par.P = 200;

%% kinetic parameters of proliferation  

par.rho_DD0 = 1.5;
% 
par.u_d = 0.05;
par.u_I = 0.1;
% 
par.G_TN4 = 8e8;
par.G_TN8 = 1e9;
% 
par.beta_TN4 = 0.5;              % 0.5 Li heliyon
par.beta_TN8 = 0.3;              % 0.3 Li heliyon
% 
par.kappa_Th = 0.6;              % 0.2  Li heliyon
par.kappa_Tr = 0.3;             % 0.05 Li heliyon
par.kappa_Tc = 0.8;              % 0.1  Li heliyon
% 
par.d_TN4 = 0.05;                 % 0.05 Li.Heliyon.2022
par.d_TN8 = 0.05;                 % 0.05 Li.Heliyon.2022
par.d_Th = 0.01;                 % 0.197 Lai.PNAS.2018,Hao.PLoSOne.2016   %%
par.d_Tr = 0.03;                 % 0.2 Lai.PNAS.2018  %%改
par.d_Tc = 0.005;                 % 0.18 Lai.PNAS.2018  %%改
% 
par.eta_Th = 1e-11;
par.eta_Tc = 5e-10;

%% kinetic parameters of tumor 
% 
par.G_X2 = 1e10;
par.G_X3 = 8e8;
par.G_X4 = 5e8;
% par.X = par.G_X1 + par.G_X2 + par.G_X3 + par.G_X4;
% 
par.beta_X1 = 0.15;
par.beta_X2 = 0.1;
par.beta_X3 = 0.1;
par.beta_X4 = 0.05;
%
par.d_X2 = 0.02;
par.d_X3 = 0.02;
par.d_X4 = 0.01;
% 
par.rho_X2X1 = 0.25;
par.rho_X3X1 = 0.25;
par.rho_X4X1 = 0.05;
par.rho_X4X2 = 0.05;
par.rho_X4X3 = 0.05;

par.Lambda_A = 2e-3;
par.Lambda_B = 2e-3;

par.K_A = 8e2;
par.K_B = 8e2;

%% kinetic parameters of pharmacokinetics 

par.A_D = 0;
par.B_D = 0;

end