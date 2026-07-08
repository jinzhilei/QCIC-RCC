function fig_6()
I = -0.1;
A = 0.1;
color_box = [0.40, 0.85, 1.00;0.60, 0.87, 0.20;1.00, 0.60, 0.00]; 
softness = 0.4;  
original = [1.00, 0.60, 0.00;0.60, 0.87, 0.20;0.40, 0.85, 1.00];
color_piont = original * (1 - softness) + softness * 1;
color_p = [1.00, 0.60, 0.00; 
             0.60, 0.87, 0.20; 
             0.40, 0.85, 1.00] *0.9;
word_size = 9;
point_size = 25;
%% load data
file='data/fig.6/PD.mat';
PD=cell2mat(struct2cell(load(file)));
file='data/fig.6/SD.mat';
SD=cell2mat(struct2cell(load(file)));
file='data/fig.6/PR.mat';
PR=cell2mat(struct2cell(load(file)));
file='data/fig.6/SR.mat';
CR=cell2mat(struct2cell(load(file)));

x=0.775+rand(1,10)*0.45;
xx=1.775+rand(1,10)*0.45;
xxx=2.775+rand(1,10)*0.45;

% plot
figure(1)
t = tiledlayout(2, 12);
nexttile([1, 3])
boxplot(CR);
hBox = findobj(gca, 'Tag', 'Box');
for i = 1:length(hBox)
    set(hBox(i), 'Color', color_box(i,:), 'LineWidth', 1.5);
    xd = get(hBox(i), 'XData');
    yd = get(hBox(i), 'YData');
    patch(xd, yd, color_box(i,:), 'FaceAlpha', 0.8, 'EdgeColor', 'none');
end
hold on
meanVals = median(CR, 1); 
xPos = 1:size(CR,2);
for i = 1:length(meanVals)
    plot(xPos(i), meanVals(i), 'o', ...
        'MarkerSize', 10, ...
        'MarkerFaceColor', color_piont(i,:), ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', 1);
end
hold on
plot(xPos, meanVals, '-', 'LineWidth', 2, 'MarkerSize', 8,'Color', [0.3, 0.3, 0.3]);
hold on
plot(x,CR(:,1),'.','Color',color_p(1,:),'MarkerSize', point_size);
hold on
plot(xx,CR(:,2),'.','Color',color_p(2,:),'MarkerSize', point_size)
hold on
plot(xxx,CR(:,3),'.','Color',color_p(3,:),'MarkerSize', point_size)
grid on

set(gca,'LineWidth',1.1);
set(gca,'FontSize',word_size);
set(gca, 'Box', 'on');
set(gca,'FontName','Airal');  % Set font name
set(gca,'XTickLabel',{'M','A','D'},'FontSize',word_size);
title('CR','FontSize',word_size);

nexttile([1, 3])
boxplot(PR);
hBox = findobj(gca, 'Tag', 'Box');
for i = 1:length(hBox)
    set(hBox(i), 'Color', color_box(i,:), 'LineWidth', 1.5);
    xd = get(hBox(i), 'XData');
    yd = get(hBox(i), 'YData');
    patch(xd, yd, color_box(i,:), 'FaceAlpha', 0.8, 'EdgeColor', 'none');
end
hold on
meanVals = median(PR, 1); 
xPos = 1 : size(PR,2);
for i = 1:length(meanVals)
    plot(xPos(i), meanVals(i), 'k-o', ...
        'MarkerSize', 10, ...
        'MarkerFaceColor', color_piont(i,:), ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', 1);
end
hold on
plot(xPos, meanVals, '-', 'LineWidth', 2, 'MarkerSize', 8,'Color',[0.3, 0.3, 0.3]);
hold on
plot(x,PR(:,1),'.','Color',color_p(1,:),'MarkerSize', point_size);
hold on
plot(xx,PR(:,2),'.','Color',color_p(2,:),'MarkerSize', point_size)
hold on
plot(xxx,PR(:,3),'.','Color',color_p(3,:),'MarkerSize', point_size)
grid on
set(gca,'LineWidth',1.1);
set(gca,'FontSize',word_size);
set(gca, 'Box', 'on');
set(gca,'FontName','Airal');  % Set font name
set(gca,'XTickLabel',{'M','A','D'},'FontSize',word_size);
title('PR','FontSize',word_size);

nexttile([1, 3])
boxplot(SD);
hBox = findobj(gca, 'Tag', 'Box');
for i = 1:length(hBox)
    set(hBox(i), 'Color', color_box(i,:), 'LineWidth', 1.5);
    xd = get(hBox(i), 'XData');
    yd = get(hBox(i), 'YData');
    patch(xd, yd, color_box(i,:), 'FaceAlpha', 0.8, 'EdgeColor', 'none');
end
hold on
meanVals = median(SD, 1); 
xPos = 1 : size(SD,2);
for i = 1:length(meanVals)
    plot(xPos(i), meanVals(i), 'k-o', ...
        'MarkerSize', 10, ...
        'MarkerFaceColor', color_piont(i,:), ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', 1);
end
hold on
plot(xPos, meanVals, '-', 'LineWidth', 2, 'MarkerSize', 8,'Color',[0.3, 0.3, 0.3]);
hold on
plot(x,SD(:,1),'.','Color',color_p(1,:),'MarkerSize', point_size);
hold on
plot(xx,SD(:,2),'.','Color',color_p(2,:),'MarkerSize', point_size)
hold on
plot(xxx,SD(:,3),'.','Color',color_p(3,:),'MarkerSize', point_size)
grid on
set(gca,'LineWidth',1.1);
set(gca,'FontSize',word_size);
set(gca, 'Box', 'on');
set(gca,'FontName','Airal');  % Set font name
set(gca,'XTickLabel',{'M','A','D'},'FontSize',word_size);
title('SD','FontSize',word_size);

nexttile([1, 3])
boxplot(PD);
hBox = findobj(gca, 'Tag', 'Box');
for i = 1:length(hBox)
    set(hBox(i), 'Color', color_box(i,:), 'LineWidth', 1.5);
    xd = get(hBox(i), 'XData');
    yd = get(hBox(i), 'YData');
    patch(xd, yd, color_box(i,:), 'FaceAlpha', 0.8, 'EdgeColor', 'none');
end
hold on
meanVals = median(PD, 1); 
xPos = 1 : size(PD,2);
for i = 1:length(meanVals)
    plot(xPos(i), meanVals(i), 'k-o', ...
        'MarkerSize', 10, ...
        'MarkerFaceColor', color_piont(i,:), ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', 1);
end
hold on
plot(xPos, meanVals, '-', 'LineWidth', 2, 'MarkerSize', 8,'Color',[0.3, 0.3, 0.3]);
hold on
plot(x,PD(:,1),'.','Color',color_box(3,:),'MarkerSize', point_size);
hold on
plot(xx,PD(:,2),'.','Color',color_box(2,:),'MarkerSize', point_size)
hold on
plot(xxx,PD(:,3),'.','Color',color_box(1,:),'MarkerSize', point_size)
grid on

set(gca,'LineWidth',1.1);
set(gca,'FontSize',word_size);
set(gca, 'Box', 'on');
set(gca,'FontName','Airal');  % Set font name
set(gca,'XTickLabel',{'M','A','D'},'FontSize',word_size);
title('PD','FontSize',word_size);
%% 
file=['data/fig.6/treat1_','Recurrence_rate_84.mat'];
Recurrence_rate_11=cell2mat(struct2cell(load(file)));

file=['data/fig.6/treat2_','Recurrence_rate_84.mat'];
Recurrence_rate_21=cell2mat(struct2cell(load(file)));

file=['data/fig.6/treat3_','Recurrence_rate_84.mat'];
Recurrence_rate_31=cell2mat(struct2cell(load(file)));

%% Processing data

Recurrence_rate_11 = flipud(Recurrence_rate_11');
Recurrence_rate_21 = flipud(Recurrence_rate_21');
Recurrence_rate_31 = flipud(Recurrence_rate_31');
Recurrence_rate_11 = round(Recurrence_rate_11, 2); 
Recurrence_rate_21 = round(Recurrence_rate_21, 2); 
Recurrence_rate_31 = round(Recurrence_rate_31, 2); 

Recurrence_rate_11_ = Recurrence_rate_21-Recurrence_rate_11;
Recurrence_rate_12_ = Recurrence_rate_31-Recurrence_rate_11;
Recurrence_rate_13_ = Recurrence_rate_31-Recurrence_rate_21;


nexttile([1, 4]) 
xvalues = {'600','1200','1800','2400'};
yvalues = {'2100','1575','1050','525'};
d = heatmap(xvalues,yvalues,Recurrence_rate_11_,'ColorbarVisible','off');
clim([I, A]);
d.Colormap = Color_RdYlBu(256);
% d.Colormap = [0 245 255; 0 229 238; 0 197 205;102 205 170;118 238 198;127 255 212;238 48 167;255 52 179]./255;
cmap = d.Colormap;
clims = d.ColorLimits;
title('Average dose VS Maximum dose');
xlabel('Atezolizumab dose');
ylabel('Bevacizumab dose');
set(gca,'FontSize',word_size);   % Set font size
set(gca,'FontName','Airal');  % Set font name

nexttile([1, 4]) 
c = heatmap(xvalues,yvalues,Recurrence_rate_12_,'ColorbarVisible','off');
c.Colormap = cmap;
c.ColorLimits = clims;
clim([I, A]);
c.Colormap = Color_RdYlBu(256);
% c.Colormap = [0 245 255; 0 229 238; 0 197 205;102 205 170;118 238 198;127 255 212;238 48 167;255 52 179]./255;
title({'Decreasing dose VS Maximum dose'});
xlabel('Atezolizumab dose');
set(gca,'FontSize',word_size);   % Set font size
set(gca,'FontName','Airal');  % Set font name

nexttile([1, 4]) 
b = heatmap(xvalues,yvalues,Recurrence_rate_13_);
b.Colormap = cmap;
b.ColorLimits = clims;
clim([I, A]);
b.Colormap = Color_RdYlBu(256);
% b.Colormap = [0 245 255; 0 229 238; 0 197 205;102 205 170;118 238 198;127 255 212;238 48 167;255 52 179]./255;
title('Decreasing dose VS Average dose');
xlabel('Atezolizumab dose');
set(gca,'FontSize',word_size);   % Set font size
set(gca,'FontName','Airal');  % Set font name
t.TileSpacing = 'compact';
t.Padding = 'compact';
set(gcf,'unit','centimeters','position',[1 1 25 13]);
print('-dpng','-r100','picture/fig.6.png')





end