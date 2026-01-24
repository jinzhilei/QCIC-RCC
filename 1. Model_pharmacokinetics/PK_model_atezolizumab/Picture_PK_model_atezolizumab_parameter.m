
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This function is used to plot the distribution histogram of parameters for the pharmacokinetics of atezolizumab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Picture_PK_model_atezolizumab_parameter()

load('./VP/VP_PK_model_atezolizumab_parameter.mat');

%% V_1_A

figure(1)

V_1_A = histogram(S(:,1));
V_1_A.NumBins = 20;
V_1_A.FaceColor = [144 238 144]./255;


%% Output image setting

% Set outer border
set(gca,'Linewidth',1);
% Set font size
set(gca,'FontSize',16);
% Set font name
set(gca,'FontName','Times New Roman');
% Set label
xlabel('L','Interpreter','latex');
ylabel('Counts','Interpreter','latex');
title({'The volume of distribution'; 'in the central compartment ($V_1$)'},'Interpreter','latex')
% Set the image name, format, and resolution
set(gcf,'unit','centimeters','position',[10 10 15 12])
print('Picture_Atezolizumab_V1','-dpng','-r600')

%% V_2_A

figure(2)

V_2_A = histogram(S(:,2));
V_2_A.NumBins = 20;
V_2_A.FaceColor = [144 238 144]./255;


%% Output image setting

% Set outer border
set(gca,'Linewidth',1);
% Set font size
set(gca,'FontSize',16);
% Set font name
set(gca,'FontName','Times New Roman');
% Set label
xlabel('L','Interpreter','latex');
ylabel('Counts','Interpreter','latex');
title({'The volume of distribution'; 'in the peripheral compartment ($V_2$)'},'Interpreter','latex')
% Set the image name, format, and resolution
set(gcf,'unit','centimeters','position',[10 10 15 12])
print('Picture_Atezolizumab_V2','-dpng','-r600')

%% Q_A

figure(3)

Q_A = histogram(S(:,3));
Q_A.NumBins = 20;
Q_A.FaceColor = [144 238 144]./255;


%% Output image setting

% Set outer border
set(gca,'Linewidth',1);
% Set font size
set(gca,'FontSize',16);
% Set font name
set(gca,'FontName','Times New Roman');
% Set label
xlabel('L/day','Interpreter','latex');
ylabel('Counts','Interpreter','latex');
title({'Capillary filtration rate between'; 'central and peripheral compartment ($Q$)'},'Interpreter','latex')
% Set the image name, format, and resolution
set(gcf,'unit','centimeters','position',[10 10 15 12])
print('Picture_Atezolizumab_Q','-dpng','-r600')


%% CI_A

figure(4)

CI_A = histogram(S(:,4));
CI_A.NumBins = 20;
CI_A.FaceColor = [144 238 144]./255;


%% Output image setting

% Set outer border
set(gca,'Linewidth',1);
% Set font size
set(gca,'FontSize',16);
% Set font name
set(gca,'FontName','Times New Roman');
% Set label
xlabel('L/day','Interpreter','latex');
ylabel('Counts','Interpreter','latex');
title({'Linear clearance rate from'; 'central compartment ($\rm{CL}$)'},'Interpreter','latex')
% Set the image name, format, and resolution
set(gcf,'unit','centimeters','position',[10 10 15 12])
print('Picture_Atezolizumab_CL','-dpng','-r600')


%% V_max_A

figure(5)

V_max_A = histogram(S(:,5));
V_max_A.NumBins = 20;
V_max_A.FaceColor = [144 238 144]./255;


%% Output image setting

% Set outer border
set(gca,'Linewidth',1);
% Set font size
set(gca,'FontSize',16);
% Set font name
set(gca,'FontName','Times New Roman');
% Set label
xlabel('mg $\cdot$ L/day','Interpreter','latex');
ylabel('Counts','Interpreter','latex');
title({'Maximum non-linear clearance rate';'($V_{max}$)'},'Interpreter','latex')
% Set the image name, format, and resolution
set(gcf,'unit','centimeters','position',[10 10 15 12])
print('Picture_Atezolizumab_V_max','-dpng','-r600')

%% K_A

figure(6)

K_A = histogram(S(:,6));
K_A.NumBins = 20;
K_A.FaceColor = [144 238 144]./255;


%% Output image setting

% Set outer border
set(gca,'Linewidth',1);
% Set font size
set(gca,'FontSize',16);
% Set font name
set(gca,'FontName','Times New Roman');
% Set label
xlabel('mg','Interpreter','latex');
ylabel('Counts','Interpreter','latex');
title({'Half-saturation of bevacizumab';'($K_M$)'},'Interpreter','latex')
% Set the image name, format, and resolution
set(gcf,'unit','centimeters','position',[10 10 15 12])
print('Picture_Atezolizumab_K_B','-dpng','-r600')