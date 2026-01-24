function picture_recurrence_rate()

%% 
Time = 360;
word_size = 13;
lengthe = 8;
hight = 15;
color = ColorMatrix();
%% Load data

file=['Recurrence_rate/treat_7/','Recurrence_rate.mat'];
Recurrence_rate_1=cell2mat(struct2cell(load(file)));
Recurrence_rate(1,:) = Recurrence_rate_1';

file=['Recurrence_rate/treat_21/','Recurrence_rate.mat'];
Recurrence_rate_2=cell2mat(struct2cell(load(file)));
Recurrence_rate(2,:) = Recurrence_rate_2';

%% Processing data

Recurrence_rate = flipud(Recurrence_rate);
Recurrence_rate = round(Recurrence_rate, 2); 
color = [cell2mat(color(1,1));cell2mat(color(1,4))];
%% Plot
% 
% figure(1)
% xvalues = {'(1200 1050)','(1800 1575)','(2400 2100)'};
% yvalues = {'21','7'};
% h = heatmap(xvalues,yvalues,Recurrence_rate);
% colormap(parula);
% title('Clinical administration');
% h.XLabel = 'Maximum dosage';
% h.YLabel = 'Medication cycle';
% set(gca,'FontSize',word_size);   % Set font size
% set(gca,'FontName','Airal');  % Set font name
% caxis([0.5, 1]);
% colorbar;
% 
% set(gcf,'unit','centimeters','position',[1 1 lengthe hight]); 
% print('-dpng','-r600','picture/treat_360.png')

figure(2)
for k = 1:2
plot(1:3,Recurrence_rate(k,:),'Color',color(k,:));
hold on
end
for k =1:2
scatter(1:3,Recurrence_rate(k,:),[],color(k,:),'filled')
hold on
for j = 1:3
    text(j-0.15,Recurrence_rate(k,j)+(-1)^(k)*0.01,[num2str(Recurrence_rate(k,j))],'FontSize', 12);
    hold on
end
end
xticks({}); % 设置x轴的刻度位置
% xticklabels({'(1200,1050)','(1800,1575)','(2400,2100)'});
legend('7 day ','21 day')
xlabel('Maximum dosing regimen','FontSize',word_size);
ylabel('Recurrence rate','FontSize',word_size);
set(gca,'XLim',[0.75 3.25]);  % Set coordinate range (X)
set(gca,'YLim',[0.55 0.8]);  % Set coordinate range (Y)
set(gca,'FontName','Airal');  % Set font name
set(gcf,'unit','centimeters','position',[1 1 lengthe hight]); 
print('-dpng','-r600','picture/pre.png')

figure(3)
drug = [1,1;1.5,1.5;2,2];
b = bar(drug);
b(1).FaceColor = cell2mat(color(4,11));
b(2).FaceColor = cell2mat(color(1,14));
xticks({}); % 设置x轴的刻度位置
yticks({}); % 设置x轴的刻度位置

set(gca,'YLim',[0 10]);  % Set coordinate range (Y)
set(gca,'FontName','Airal');  % Set font name
set(gcf,'unit','centimeters','position',[1 1 lengthe hight]); 
print('-dpng','-r600','picture/pre_drug.png')







end