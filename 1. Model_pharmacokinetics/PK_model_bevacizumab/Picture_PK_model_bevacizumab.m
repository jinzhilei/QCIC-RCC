
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This function is used to plot the dynamic change curve of bevacizumab concentration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Picture_PK_model_bevacizumab()

sample = 1000;

load('./VP/VP_central.mat');

load('./VP/VP_peripheral.mat');

%% Plot the solution curve for each PK model

% for i=1:sample
%     P = semilogy(VP_central(:,1),VP_central(:,i+1));
%     P.LineWidth = 3; 
%     P.Color = [30 144 255]./255;
%     hold on;
% end

%% Plot the 95% confidence interval and average curve of the PK model

% 100% Confidence interval
% P_up = semilogy(VP_central(:,1),max(VP_central(:,2:sample+1),[],2),':');
% P_up.LineWidth = 2; 
% P_up.Color = [135 206 255]./255;
% hold on

% P_down = semilogy(VP_central(:,1),min(VP_central(:,2:sample+1),[],2),':');
% P_down.LineWidth = 2; 
% P_down.Color = [135 206 255]./255;
% hold on

figure(1)

S_95CI = sort(VP_central(:,2:sample+1),2);

P_up_central = semilogy(VP_central(:,1),S_95CI(:,26),':');
P_up_central.LineWidth = 2; 
P_up_central.Color = [135 206 255]./255;
hold on

P_down_central = semilogy(VP_central(:,1),S_95CI(:,975),':');
P_down_central.LineWidth = 2; 
P_down_central.Color = [135 206 255]./255;
hold on


t_line_central=[VP_central(:,1);flipud(VP_central(:,1))];
drug_line_central=[S_95CI(:,26);flipud(S_95CI(:,975))];
P_fill_central=fill(t_line_central, drug_line_central,[135 206 255]./255,'FaceAlpha',0.25,'EdgeColor','none');
hold on

P_mean_central = semilogy(VP_central(:,1),mean(VP_central(:,2:sample+1),2));
P_mean_central.LineWidth = 2; 
P_mean_central.Color = [30 144 255]./255;
grid on;

%% Plot experimental data points

x_errorbar = [4 6 12 20 26 30 35 43 47 50 55 65];
y_errorbar = [265 185 95 50 350 240 150 75 400 300 190 70];
yneg = [135 115 55 30 150 90 70 45 200 150 110 50];
ypos = [135 115 55 30 150 90 70 45 200 150 110 50];
P_errorbar = errorbar(x_errorbar,y_errorbar,yneg,ypos,'o','Color',[255 105 180]./255,'MarkerSize',8,'MarkerEdgeColor',[255 105 180]./255,'MarkerFaceColor',[255 105 180]./255);

%% Output image setting

% Set coordinate range
set(gca,'XLim',[1 360]);
set(gca,'YLim',[10 1400]);
% Set outer border
set(gca,'Linewidth',1);
% Set font size
set(gca,'FontSize',14);
% Set font name
set(gca,'FontName','Times New Roman');
% Set legend
legend([P_fill_central P_mean_central P_errorbar],{' 95% Confidence interval',' Mean (numerical calculation)',' Clinically measured mean ± SD'})
% Set label
xlabel('Time (day)');
ylabel('Bevacizumab concentration ($\mu$g/mL)','Interpreter','latex');
% Set the image name, format, and resolution
set(gcf,'unit','centimeters','position',[2 2 50 10])
print('Picture_Bevacizumab_serum_concentration','-dpng','-r600')

figure(2)

S_95CI_peripheral = sort(VP_peripheral(:,2:sample+1),2);

P_up_peripheral = semilogy(VP_peripheral(:,1),S_95CI_peripheral(:,26),':');
P_up_peripheral.LineWidth = 2; 
P_up_peripheral.Color = [135 206 255]./255;
hold on

P_down_peripheral = semilogy(VP_peripheral(:,1),S_95CI_peripheral(:,975),':');
P_down_peripheral.LineWidth = 2; 
P_down_peripheral.Color = [135 206 255]./255;
hold on


t_line_peripheral=[VP_peripheral(:,1);flipud(VP_peripheral(:,1))];
drug_line_peripheral=[S_95CI_peripheral(:,26);flipud(S_95CI_peripheral(:,975))];
P_fill_peripheral=fill(t_line_peripheral, drug_line_peripheral,[135 206 255]./255,'FaceAlpha',0.25,'EdgeColor','none');
hold on

P_mean_peripheral = semilogy(VP_peripheral(:,1),mean(VP_peripheral(:,2:sample+1),2));
P_mean_peripheral.LineWidth = 2; 
P_mean_peripheral.Color = [30 144 255]./255;
grid on;

%% Output image setting

% Set coordinate range
set(gca,'XLim',[1 360]);
set(gca,'YLim',[1 600]);
% Set outer border
set(gca,'Linewidth',1);
% Set font size
set(gca,'FontSize',14);
% Set font name
set(gca,'FontName','Times New Roman');
% Set legend
legend([P_fill_peripheral P_mean_peripheral],{' 95% Confidence interval',' Mean (numerical calculation)'})
% Set label
xlabel('Time (day)');
ylabel('Bevacizumab peripheral concentration ($\mu$g/mL)','Interpreter','latex');
% Set the image name, format, and resolution
set(gcf,'unit','centimeters','position',[10 10 50 12])
print('Picture_Bevacizumab_peripheral_concentration','-dpng','-r600')

