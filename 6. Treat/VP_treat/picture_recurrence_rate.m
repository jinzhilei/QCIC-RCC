function picture_recurrence_rate()

%% 
Time = 360;
word_size = 12;
lengthe = 50;
hight = 30;
xlim = 1;
%% Load data

file=['Recurrence_rate/treat',num2str(1),'_','Recurrence_rate_',num2str(Time),'.mat'];
Recurrence_rate_1=cell2mat(struct2cell(load(file)));

file=['Recurrence_rate/treat',num2str(2),'_','Recurrence_rate_',num2str(Time),'.mat'];
Recurrence_rate_2=cell2mat(struct2cell(load(file)));

file=['Recurrence_rate/treat',num2str(3),'_','Recurrence_rate_',num2str(Time),'.mat'];
Recurrence_rate_3=cell2mat(struct2cell(load(file)));

%% Processing data

Recurrence_rate_1 = flipud(Recurrence_rate_1');
Recurrence_rate_2 = flipud(Recurrence_rate_2');
Recurrence_rate_3 = flipud(Recurrence_rate_3');
Recurrence_rate_1 = round(Recurrence_rate_1, 2); 
Recurrence_rate_2 = round(Recurrence_rate_2, 2); 
Recurrence_rate_3 = round(Recurrence_rate_3, 2); 

%% Plot

figure(1)
tiledlayout(2,3);
nexttile
xvalues = {'300','600','900','1200','1500','1800','2100','2400'};
yvalues = {'2100','1837','1575','1312','1050','787','525','262'};
d = heatmap(xvalues,yvalues,Recurrence_rate_1,'ColorbarVisible','off');
colormap(parula);
cmap = d.Colormap;
clims = d.ColorLimits;
title('Clinical administration');
xlabel('Atezolizumab dose');
ylabel('Bevacizumab dose');
set(gca,'FontSize',word_size);   % Set font size
set(gca,'FontName','Airal');  % Set font name


nexttile
xvalues = {'300','600','900','1200','1500','1800','2100','2400'};
yvalues = {'2100','1837','1575','1312','1050','787','525','262'};
c = heatmap(xvalues,yvalues,Recurrence_rate_2,'ColorbarVisible','off');
c.Colormap = cmap;
c.ColorLimits = clims;
title('Sustained dosing');
xlabel('Atezolizumab dose');
ylabel('Bevacizumab dose');
set(gca,'FontSize',word_size);   % Set font size
set(gca,'FontName','Airal');  % Set font name

nexttile
xvalues = {'300','600','900','1200','1500','1800','2100','2400'};
yvalues = {'2100','1837','1575','1312','1050','787','525','262'};
b = heatmap(xvalues,yvalues,Recurrence_rate_3);
b.Colormap = cmap;
b.ColorLimits = clims;

clim([0.4, 1]);
colormap(parula);
xlabel('Atezolizumab dose');
ylabel('Bevacizumab dose');
title('Waning treatment');
set(gca,'FontSize',word_size);   % Set font size
set(gca,'FontName','Airal');  % Set font name

nexttile
[C,p]=contour(flipud(Recurrence_rate_1),'LineWidth',2); 
clabel(C,p,'FontSize',word_size); 
hold on
plot(4,4,'r*','MarkerSize', 15, 'LineWidth', 3)
xlabel('Atezolizumab dose');
ylabel('Bevacizumab dose');
title('Clinical administration');

set(gca,'FontSize',word_size);   % Set font size
set(gca,'FontName','Airal');  % Set font name
set(gca, 'XTickLabel', xvalues);
set(gca, 'YTickLabel', sort([2100,1837,1575,1312,1050,787,525,262]));
set(gca,'ZLim',[0 xlim]);  


nexttile
[C,p]=contour(flipud(Recurrence_rate_2),'LineWidth',2); 
clabel(C,p,'FontSize',word_size); 
hold on
plot(4,4,'r*','MarkerSize', 15, 'LineWidth', 3)
xlabel('Atezolizumab dose');
ylabel('Bevacizumab dose');
title('Sustained dosing');
set(gca,'FontSize',word_size);   % Set font size
set(gca,'FontName','Airal');  % Set font name
set(gca, 'XTickLabel', xvalues);
set(gca, 'YTickLabel', sort([2100,1837,1575,1312,1050,787,525,262]));
set(gca,'ZLim',[0 xlim]);  % Set coordinate range (X)

nexttile
[C,p]=contour(flipud(Recurrence_rate_3),'LineWidth',2); 
clabel(C,p,'FontSize',word_size); 
hold on
plot(4,4,'r*','MarkerSize', 15, 'LineWidth', 3)
xlabel('Atezolizumab dose');
ylabel('Bevacizumab dose');
title('Waning treatment');
set(gca,'FontSize',word_size);   % Set font size
set(gca,'FontName','Airal');  % Set font name
set(gca, 'XTickLabel', xvalues);
set(gca, 'YTickLabel', sort([2100,1837,1575,1312,1050,787,525,262]));
set(gca,'ZLim',[0 xlim]);  % Set coordinate range (X)
colorbar
% 
set(gcf,'unit','centimeters','position',[1 1 lengthe hight]); 
print('-dpng','-r600','picture/recurrence_all_time_1.png')


end