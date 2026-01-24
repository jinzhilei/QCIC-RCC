
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This function is used to plot the dynamic change curve of atezolizumab concentration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Picture_PK_model_atezolizumab()

sample = 1000;

load('./VP/VP_central.mat');

load('./VP/VP_peripheral.mat');

%% Plot the solution curve for each PK model

% for i=1:sample
%     P = semilogy(VP_central(:,1),VP_central(:,i+1));
%     P.LineWidth = 3; 
%     P.Color = [144 238 144]./255;
%     hold on;
% end

%% Plot the 95% confidence interval and average curve of the PK model

% 100% Confidence interval
% P_up = semilogy(VP_central(:,1),max(VP_central(:,2:sample+1),[],2),':');
% P_up.LineWidth = 2; 
% P_up.Color = [144 238 144]./255;
% hold on

% P_down = semilogy(VP_central(:,1),min(VP_central(:,2:sample+1),[],2),':');
% P_down.LineWidth = 2; 
% P_down.Color = [144 238 144]./255;
% hold on

figure(1)

[S_95CI_central,I] = sort(VP_central(:,2:sample+1),2);

P_up_central = semilogy(VP_central(:,1),S_95CI_central(:,26),':');
P_up_central.LineWidth = 2; 
P_up_central.Color = [144 238 144]./255;
hold on

P_down_central = semilogy(VP_central(:,1),S_95CI_central(:,975),':');
P_down_central.LineWidth = 2; 
P_down_central.Color = [144 238 144]./255;
hold on


t_line_central=[VP_central(:,1);flipud(VP_central(:,1))];
drug_line_central=[S_95CI_central(:,26);flipud(S_95CI_central(:,975))];
P_fill_central=fill(t_line_central, drug_line_central,[144 238 144]./255,'FaceAlpha',0.25,'EdgeColor','none');
hold on

P_mean_central = semilogy(VP_central(:,1),mean(VP_central(:,2:sample+1),2));
P_mean_central.LineWidth = 2; 
P_mean_central.Color = [50 205 50]./255;
grid on;

%% Plot experimental data points

x_errorbar = [3 7 14 20];
y_errorbar = [325 150 120 90];
yneg = [75 50 45 30];
ypos = [75 50 45 30];
P_errorbar = errorbar(x_errorbar,y_errorbar,yneg,ypos,'o','Color',[255 165 0]./255,'MarkerSize',8,'MarkerEdgeColor',[255 165 0]./255,'MarkerFaceColor',[255 165 0]./255);

%% Output image setting

% Set coordinate range
set(gca,'XLim',[1 150]);
set(gca,'YLim',[10 2000]);
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
ylabel('Atezolizumab concentration ($\mu$g/mL)','Interpreter','latex');
% Set the image name, format, and resolution
set(gcf,'unit','centimeters','position',[2 2 50 10])
print('Picture_Atezolizumab_serum_concentration','-dpng','-r600')

figure(2)

S_95CI_peripheral = sort(VP_peripheral(:,2:sample+1),2);

P_up_peripheral = semilogy(VP_peripheral(:,1),S_95CI_peripheral(:,26),':');
P_up_peripheral.LineWidth = 2; 
P_up_peripheral.Color = [144 238 144]./255;
hold on

P_down_peripheral = semilogy(VP_peripheral(:,1),S_95CI_peripheral(:,975),':');
P_down_peripheral.LineWidth = 2; 
P_down_peripheral.Color = [144 238 144]./255;
hold on


t_line_peripheral=[VP_peripheral(:,1);flipud(VP_peripheral(:,1))];
drug_line_peripheral=[S_95CI_peripheral(:,26);flipud(S_95CI_peripheral(:,975))];
P_fill_peripheral=fill(t_line_peripheral, drug_line_peripheral,[144 238 144]./255,'FaceAlpha',0.25,'EdgeColor','none');
hold on

P_mean_peripheral = semilogy(VP_peripheral(:,1),mean(VP_peripheral(:,2:sample+1),2));
P_mean_peripheral.LineWidth = 2; 
P_mean_peripheral.Color = [50 205 50]./255;
grid on;

%% Output image setting

% Set coordinate range
set(gca,'XLim',[1 360]);
set(gca,'YLim',[1 600]);
% Set outer border
set(gca,'Linewidth',1);
% Set font size
set(gca,'FontSize',12);
% Set font name
set(gca,'FontName','Times New Roman');
% Set legend
legend([P_fill_peripheral P_mean_peripheral],{' 95% Confidence interval',' Mean (numerical calculation)'})
% Set label
xlabel('Time (day)');
ylabel('Atezolizumab peripheral concentration ($\mu$g/mL)','Interpreter','latex');
% Set the image name, format, and resolution
set(gcf,'unit','centimeters','position',[2 2 30 12])
print('Picture_Atezolizumab_peripheral_concentration','-dpng','-r600')

%% Save load

file='./VP/VP_central_sort.mat';
save(file,'I');
