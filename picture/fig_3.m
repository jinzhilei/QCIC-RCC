function fig_3()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Plot: This file is used to plot figure 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sample = 1000;


VP_central_A = cell2mat(struct2cell(load('./data/fig.3/A/VP_central_A.mat')));
S_A = cell2mat(struct2cell(load('./data/fig.3/A/VP_PK_model_atezolizumab_parameter.mat')));
VP_central_B = cell2mat(struct2cell(load('./data/fig.3/B/VP_central_B.mat')));
S_B = cell2mat(struct2cell(load('./data/fig.3/B/VP_PK_model_bevacizumab_parameter.mat')));
%% Plot the 95% confidence interval and average curve of the PK model

figure(1)
subplot(2,2,1)

%% A
[S_95CI_central,~] = sort(VP_central_A(:,2:sample+1),2);

P_up_central = semilogy(VP_central_A(:,1),S_95CI_central(:,26),':');
P_up_central.LineWidth = 2; 
P_up_central.Color = [144 238 144]./255;
hold on

P_down_central = semilogy(VP_central_A(:,1),S_95CI_central(:,975),':');
P_down_central.LineWidth = 2; 
P_down_central.Color = [144 238 144]./255;
hold on


t_line_central=[VP_central_A(:,1);flipud(VP_central_A(:,1))];
drug_line_central=[S_95CI_central(:,26);flipud(S_95CI_central(:,975))];
P_fill_central=fill(t_line_central, drug_line_central,[144 238 144]./255,'FaceAlpha',0.25,'EdgeColor','none');
hold on

P_mean_central = semilogy(VP_central_A(:,1),mean(VP_central_A(:,2:sample+1),2));
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
set(gca,'XLim',[1 64]);
set(gca,'YLim',[10 10000]);
% Set outer border
set(gca,'Linewidth',1);
% Set font size
set(gca,'FontSize',10);
% Set legend
legend([P_fill_central P_mean_central P_errorbar],{' 95% Confidence interval',' Mean (numerical calculation)',' Clinically measured mean ± SD'},'FontSize',10)
legend boxoff
% Set label
xlabel('Time (day)');
ylabel('Atezolizumab (\mug/mL)');
set(gca,'FontName','Arial');

%% V_1_A
subplot(4,6,4)
V_1_A = histogram(S_A(:,1));
V_1_A.NumBins = 20;
V_1_A.FaceColor = [144 238 144]./255;

%% Output image setting

% Set outer border
set(gca,'Linewidth',1);
% Set font size
set(gca,'FontSize',10);
% Set label
xlabel('V_1 (L)');
ylabel('Counts');
% Set font name
set(gca,'FontName','Arial');
grid on


%% V_2_A

subplot(4,6,5)
V_2_A = histogram(S_A(:,2));
V_2_A.NumBins = 20;
V_2_A.FaceColor = [144 238 144]./255;


%% Output image setting

% Set outer border
set(gca,'Linewidth',1);
% Set font size
set(gca,'FontSize',10);
% Set label
xlabel('V_2 (L)');
% Set font name
set(gca,'FontName','Arial');
grid on

%% Q_A

subplot(4,6,6)

Q_A = histogram(S_A(:,3));
Q_A.NumBins = 20;
Q_A.FaceColor = [144 238 144]./255;


%% Output image setting

% Set outer border
set(gca,'Linewidth',1);
% Set font size
set(gca,'FontSize',10);
% Set label
xlabel('Q (L/day)');
% Set font name
set(gca,'FontName','Arial');
grid on

%% CI_A

subplot(4,6,10)

CI_A = histogram(S_A(:,4));
CI_A.NumBins = 20;
CI_A.FaceColor = [144 238 144]./255;


%% Output image setting

% Set outer border
set(gca,'Linewidth',1);
% Set font size
set(gca,'FontSize',10);
% Set label
xlabel('\rm{CL} (L/day)');
ylabel('Counts');
% Set font name
set(gca,'FontName','Arial');
grid on

%% V_max_A

subplot(4,6,11)

V_max_A = histogram(S_A(:,5));
V_max_A.NumBins = 20;
V_max_A.FaceColor = [144 238 144]./255;


%% Output image setting

% Set outer border
set(gca,'Linewidth',1);
% Set font size
set(gca,'FontSize',10);
% Set label
xlabel(' V_{max} (mg · L/day)');
% Set font name
set(gca,'FontName','Arial');
grid on

%% K_A

subplot(4,6,12)

K_A = histogram(S_A(:,6));
K_A.NumBins = 20;
K_A.FaceColor = [144 238 144]./255;


%% Output image setting

% Set outer border
set(gca,'Linewidth',1);
% Set font size
set(gca,'FontSize',10);
% Set label
xlabel('K_M (mg)');
% Set font name
set(gca,'FontName','Arial');
grid on

%% B
subplot(2,2,3)
S_95CI = sort(VP_central_B(:,2:sample+1),2);

P_up_central = semilogy(VP_central_B(:,1),S_95CI(:,26),':');
P_up_central.LineWidth = 2; 
P_up_central.Color = [135 206 255]./255;
hold on

P_down_central = semilogy(VP_central_B(:,1),S_95CI(:,975),':');
P_down_central.LineWidth = 2; 
P_down_central.Color = [135 206 255]./255;
hold on


t_line_central=[VP_central_B(:,1);flipud(VP_central_B(:,1))];
drug_line_central=[S_95CI(:,26);flipud(S_95CI(:,975))];
P_fill_central=fill(t_line_central, drug_line_central,[135 206 255]./255,'FaceAlpha',0.25,'EdgeColor','none');
hold on

P_mean_central = semilogy(VP_central_B(:,1),mean(VP_central_B(:,2:sample+1),2));
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
set(gca,'XLim',[1 64]);
set(gca,'YLim',[10 10000]);
% Set outer border
set(gca,'Linewidth',1);
% Set font size
set(gca,'FontSize',10);
% Set font name
% Set legend
legend([P_fill_central P_mean_central P_errorbar],{' 95% Confidence interval',' Mean (numerical calculation)',' Clinically measured mean ± SD'},'FontSize',10)
legend boxoff
% Set label
xlabel('Time (day)');
ylabel('Bevacizumab (\mu g/mL)');
set(gca,'FontName','Arial');


%% V_1_B

subplot(4,6,16)
V_1_B = histogram(S_B(:,1));
V_1_B.NumBins = 20;
V_1_B.FaceColor = [135 206 255]./255;


%% Output image setting

% Set outer border
set(gca,'Linewidth',1);
% Set font size
set(gca,'FontSize',10);
% Set label
xlabel('V_1 (L)');
ylabel('Counts');
set(gca,'FontName','Arial');
grid on

%% V_2_B

subplot(4,6,17)
V_2_B = histogram(S_B(:,2));
V_2_B.NumBins = 20;
V_2_B.FaceColor = [135 206 255]./255;


%% Output image setting

% Set outer border
set(gca,'Linewidth',1);
% Set font size
set(gca,'FontSize',10);
% Set label
xlabel('V_2 (L)');
set(gca,'FontName','Arial');
grid on

%% Q_B

subplot(4,6,18)

Q_B = histogram(S_B(:,3));
Q_B.NumBins = 20;
Q_B.FaceColor = [135 206 255]./255;


%% Output image setting

% Set outer border
set(gca,'Linewidth',1);
% Set font size
set(gca,'FontSize',10);
% Set label
xlabel('Q (L/day)');
set(gca,'FontName','Arial');
grid on


%% CI_B

subplot(4,6,22)

CI_B = histogram(S_B(:,4));
CI_B.NumBins = 20;
CI_B.FaceColor = [135 206 255]./255;


%% Output image setting

% Set outer border
set(gca,'Linewidth',1);
% Set font size
set(gca,'FontSize',10);
% Set label
xlabel('CI (L/day)');
ylabel('Counts');
set(gca,'FontName','Arial');
grid on


%% V_max_B

subplot(4,6,23)

V_max_B = histogram(S_B(:,5));
V_max_B.NumBins = 20;
V_max_B.FaceColor = [135 206 255]./255;


%% Output image setting

% Set outer border
set(gca,'Linewidth',1);
% Set font size
set(gca,'FontSize',10);
% Set label
xlabel('V_{max} (mg \cdot L/day)');
set(gca,'FontName','Arial');
grid on

%% K_B

subplot(4,6,24)

K_B = histogram(S_B(:,6));
K_B.NumBins = 20;
K_B.FaceColor = [135 206 255]./255;


%% Output image setting

% Set outer border
set(gca,'Linewidth',1);
% Set font size
set(gca,'FontSize',10);
% Set label
xlabel('K_M (mg)');
set(gca,'FontName','Arial');
grid on 

% Set the image name, format, and resolution
set(gcf,'unit','centimeters','position',[2 2 40 20])

print('picture/fig.3.png','-dpng','-r100')


