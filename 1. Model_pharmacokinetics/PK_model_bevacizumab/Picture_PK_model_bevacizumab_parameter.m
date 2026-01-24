
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This function is used to plot the distribution histogram of parameters for the pharmacokinetics of bevacizumab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Picture_PK_model_bevacizumab_parameter()

load('./VP/VP_PK_model_bevacizumab_parameter.mat');

%% V_1_B

figure(1)

V_1_B = histogram(S(:,1));
V_1_B.NumBins = 20;
V_1_B.FaceColor = [135 206 255]./255;


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
print('Picture_Bevacizumab_V1','-dpng','-r600')

%% V_2_B

figure(2)

V_2_B = histogram(S(:,2));
V_2_B.NumBins = 20;
V_2_B.FaceColor = [135 206 255]./255;


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
print('Picture_Bevacizumab_V2','-dpng','-r600')

%% Q_B

figure(3)

Q_B = histogram(S(:,3));
Q_B.NumBins = 20;
Q_B.FaceColor = [135 206 255]./255;


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
print('Picture_Bevacizumab_Q','-dpng','-r600')


%% CI_B

figure(4)

CI_B = histogram(S(:,4));
CI_B.NumBins = 20;
CI_B.FaceColor = [135 206 255]./255;


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
print('Picture_Bevacizumab_CL','-dpng','-r600')


%% V_max_B

figure(5)

V_max_B = histogram(S(:,5));
V_max_B.NumBins = 20;
V_max_B.FaceColor = [135 206 255]./255;


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
print('Picture_Bevacizumab_V_max','-dpng','-r600')

%% K_B

figure(6)

K_B = histogram(S(:,6));
K_B.NumBins = 20;
K_B.FaceColor = [135 206 255]./255;


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
print('Picture_Bevacizumab_K_B','-dpng','-r600')