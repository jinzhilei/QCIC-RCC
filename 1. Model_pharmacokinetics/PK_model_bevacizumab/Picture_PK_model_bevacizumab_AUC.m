
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This function is used to plot the AUC curve

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Picture_PK_model_bevacizumab_AUC()

sample = 1000;

load('./VP/VP_central.mat');

AUC_individual = zeros(21,sample);

for i=1:sample

t = VP_central(1:22,1);
C = VP_central(1:22,i+1);

AUC_individual(:,i) = Numerical_integration(t,C);

% plot(1:1:21,AUC_individual(:,i));
% hold on

end

AUC_95CI = sort(AUC_individual,2);

figure(1)

P_up = plot(1:1:21,AUC_95CI(:,975),':');
P_up.LineWidth = 2; 
P_up.Color = [135 206 255]./255;
hold on

P_down = plot(1:1:21,AUC_95CI(:,26),':');
P_down.LineWidth = 2; 
P_down.Color = [135 206 255]./255;
hold on

t_line=[1:1:21 fliplr(1:1:21)];
drug_line=[AUC_95CI(:,26);flipud(AUC_95CI(:,975))];
P_fill=fill(t_line, drug_line,[135 206 255]./255,'FaceAlpha',0.25,'EdgeColor','none');
hold on

P_mean = plot(1:1:21,mean(AUC_individual,2));
P_mean.LineWidth = 2; 
P_mean.Color = [30 144 255]./255;

%% Output image setting

% Set coordinate range
set(gca,'XLim',[0 22]);
set(gca,'YLim',[0 7000]);
% Set outer border
set(gca,'Linewidth',1);
% Set font size
set(gca,'FontSize',14);
% Set font name
set(gca,'FontName','Times New Roman');
% Set legend
legend([P_fill P_mean],{' 95% Confidence interval',' Mean (numerical calculation)'},'Location','Northwest')
% Set label
xlabel('Time (day)');
ylabel('AUC');
title('Numerical result type');
% Set the image name, format, and resolution
set(gcf,'unit','centimeters','position',[10 10 15 12])
print('Picture_Bevacizumab_AUC_Cycle1','-dpng','-r600')


figure(2)

AUC_time = 1:1:21;
AUC_mean = mean(AUC_individual,2);
AUC_up = AUC_95CI(:,975);
AUC_down = AUC_95CI(:,26);

% cftool(AUC_time,AUC_mean)
AUC_function_mean = 1108.7*log(AUC_time)+879.8939;

% cftool(AUC_time,AUC_up)
AUC_function_up = 1586.7*log(AUC_time)+632.5774;

% cftool(AUC_time,AUC_down)
AUC_function_down = 680.1793*log(AUC_time)+1038.7;

P_1 = plot(AUC_time,AUC_function_up,'-','LineWidth',2,'Color',[0 61 245]./255);
hold on
P_2 = plot(AUC_time,AUC_function_mean,'-','LineWidth',2,'Color',[30 144 255]./255);
hold on
P_3 = plot(AUC_time,AUC_function_down,'-','LineWidth',2,'Color',[112 219 255]./255);
hold on

for i=1:21
    P_4 = plot(i,AUC_up(i),'o','Color',[0 61 245]./255);
    hold on
    P_5 = plot(i,AUC_mean(i),'o','Color',[30 144 255]./255);
    hold on
    P_6 = plot(i,AUC_down(i),'o','Color',[112 219 255]./255);
    hold on
end

%% Output image setting

% Set coordinate range
set(gca,'XLim',[0 22]);
set(gca,'YLim',[0 7000]);
% Set outer border
set(gca,'Linewidth',1);
% Set font size
set(gca,'FontSize',14);
% Set font name
set(gca,'FontName','Times New Roman');
% Set legend
legend([P_1 P_2 P_3 P_4 P_5 P_6],{'97.5th percentile','Mean','2.5th percentile','97.5th data','mean data','2.5th data'},'Location','Northwest')
% Set label
xlabel('Time (day)');
ylabel('AUC');
title('Function result type');
% Set the image name, format, and resolution
set(gcf,'unit','centimeters','position',[10 10 15 12])
print('Picture_Bevacizumab_AUC_Log','-dpng','-r600')

end
