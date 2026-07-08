function fig_8()
n=100;
t=360;
TT=1:1:t;
T=TT./30;
I = 70;
A = 100;
day = [0,90,180,270,360];
F_m = zeros(1,5);
F_n = [0,4.8,5,6.7,3.4];
file_response='./data/fig.8/OS_ave.mat';
Q_ave = cell2mat(struct2cell(load(file_response)));
for i = 2:5
    if day(i-1) == 0
        F_m(i) = Q_ave(day(i-1)+1)-Q_ave(day(i));
    else
        F_m(i) = Q_ave(day(i-1))-Q_ave(day(i));
    end
end
Title = [{'Tc'},{'Tr/Tc'},{'Tc/TAM'},{'Th/Tc'}];
OS_ave = cell(8, 1);
OS = cell(8,1);
for i =1:8
    file = ['./data/fig.8/OS_ave',num2str(i),'.mat'];
    OS_ave{i,1}= cell2mat(struct2cell(load(file)));
    file = ['./data/fig.8/OS',num2str(i),'.mat'];
    OS{i,1}= cell2mat(struct2cell(load(file)));
end
nSamples = size(OS_ave{1,1}, 2);
file_response='./data/fig.8/OS.mat';
Q = cell2mat(struct2cell(load(file_response)));
file_response='./data/fig.8/OS_ave.mat';
Q_ave = cell2mat(struct2cell(load(file_response)));

subplot(2,3,1)

plot(T,Q_ave,'-','color',[0 205 102]./255,'linewidth',3);
hold on
%% Clin Data
B1 = load('data/fig.8/1.Rini.2019.dat');
%% curve plotting
plot(B1(:,1),B1(:,2),'-','color',[30 144 255]./255,'linewidth',3);
grid on
%% Output image setting
% Set coordinate range
set(gca,'XLim',[0 12]);
set(gca,'YLim',[I A]);
% Set font size
set(gca,'FontSize',16);
% Set font name
set(gca,'FontName','Arial');
% Set label
xlabel('Time (months)');
ylabel('OS (%)');
% Set legend
legend('Model simulation','Rini et al. 2019','Location','southeast');%'95% Confidence Interval',
legend('boxoff')

subplot(2,3,2)
plot(day./30,F_m,'-','color',[0 205 102]./255,'linewidth',3);
hold on
plot(day./30,F_n,'-','color',[30 144 255]./255,'linewidth',3);
hold on
plot(day./30,F_m,'o','MarkerEdgeColor',[1,1,1],'MarkerFaceColor',[0 205 102]./255,'MarkerSize',12,'LineWidth',1.25,'DisplayName','PD');
hold on
plot(day./30,F_n,'o','MarkerEdgeColor',[1,1,1],'MarkerFaceColor',[30 144 255]./255,'MarkerSize',12,'LineWidth',1.25,'DisplayName','SD');
grid on;
set(gca,'XLim',[0 12]);
set(gca,'YLim',[0 10]);
% Set font size
set(gca,'FontSize',16);
% Set font name
set(gca,'FontName','Arial');
% Set label
xlabel('Time (months)');
ylabel('Frequency')

for i =1:4
    mean_data_1 = mean(OS{2*i-1,1}, 2);   % 每行的均值，得到31×1的向量
    std_data_1 = std(OS{2*i-1}, 0, 2);  % 每行的标准差，31×1
    sem = std_data_1 / sqrt(nSamples);  % 标准误（standard error of mean）
    alpha = 0.05;                     % 显著性水平
    t_crit = tinv(1 - alpha/2, nSamples - 1);  % t 临界值（双尾）
    CI_1_lower = mean_data_1 - t_crit * sem;
    CI_1_upper = mean_data_1 + t_crit * sem;
    mean_data_2 = mean(OS{2*i,1}, 2);   % 每行的均值，得到31×1的向量
    std_data_2 = std(OS{2*i,1}, 0, 2);  % 每行的标准差，31×1
    sem = std_data_2 / sqrt(nSamples);  % 标准误（standard error of mean）
    alpha = 0.05;                     % 显著性水平
    t_crit = tinv(1 - alpha/2, nSamples - 1);  % t 临界值（双尾）
    CI_2_lower = mean_data_2 - t_crit * sem;
    CI_2_upper = mean_data_2 + t_crit * sem;

    subplot(2,3,i+2)

    fill([T, fliplr(T)], ...
        [CI_1_upper', fliplr(CI_1_lower')], ...
        [0 255 127]./255, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    hold on
    fill([T, fliplr(T)], ...
        [CI_2_upper', fliplr(CI_2_lower')], ...
        [0 255 127]./255, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    hold on
    h1 = plot(T,mean_data_1,'-','Color',[0 205 102]./255,'LineWidth',3);
    hold on
    h2 = plot(T,mean_data_2,':','Color',[0 205 102]./255,'LineWidth',3);
    grid on;

    %% Output image setting
    % Set coordinate range
    set(gca,'XLim',[0 12]);
    set(gca,'YLim',[I A]);
    % Set font size
    set(gca,'FontSize',17);
    % Set font name
    set(gca,'FontName','Arial');
    % Set label
    xlabel('Time (months)');
    ylabel('OS (%)');
    title(Title(i));
    % Set legend
    legend([h1, h2], {'high','low'}' ,'Location','southwest');%'95% Confidence Interval',
    legend('boxoff')
end

set(gcf,'unit','centimeters','position',[3 4 35 20])
orient landscape;
exportgraphics(gcf, 'picture/OS.pdf', 'ContentType', 'vector');

