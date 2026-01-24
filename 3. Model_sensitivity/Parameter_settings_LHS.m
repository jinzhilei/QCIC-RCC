theta_D0 = 1e8;             
theta_TN4 = 5e8;           
theta_TN8 = 5e8;        
theta_M = 1e8; 
kappa_Th = 0.6;              
kappa_Tr = 0.55;           
kappa_Tc = 0.8;            
beta_Th = 0.4;              
beta_Tr = 0.5;             
beta_Tc = 0.43;  


% Parameter Labels

PRCC_var={'$$\theta_{D_0}$$' '$$\theta_{T_{N4}}$$' '$$\theta_{T_{N8}}$$' '$$\theta_M$$'...
     '$$\kappa_{T_h}$$' '$$\kappa_{T_r}$$' '$$\kappa_{T_c}$$' '$$\beta_{T_h}$$'...
     '$$\beta_{T_r}$$' '$$\beta_{T_c}$$'};% 

%% TIME SPAN OF THE SIMULATION

t_end=300; % length of the simulations
tspan=(0:1:t_end);   % time points where the output is calculated
time_points=1:300; % time points of interest for the US analysis

% Initial condition for the ODE model
 
y0=initialzation_basic;
