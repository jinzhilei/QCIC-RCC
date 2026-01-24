function picture_bar_error() 
color = ColorMatrix();
color_size = [7,10,4,1];
file='./VP_A_bar/CPSP';
CPSP=cell2mat(struct2cell(load(file)));
CPSP = sort(CPSP);
CR = mean(CPSP(:,4));
PR = mean(CPSP(:,3));
SD = mean(CPSP(:,2));
PD = mean(CPSP(:,1));

CR_real = 6.06;
PR_real = 19.19;
SD_real = 36.36;
PD_real = 38.38;

CR_error_max = CPSP(975,4)-CR;
CR_error_min = CR-CPSP(1,4);
PR_error_max = CPSP(975,3)-PR;
PR_error_min = PR-CPSP(1,3);
SD_error_max = CPSP(975,2)-SD;
SD_error_min = SD-CPSP(1,2);
PD_error_max = CPSP(975,1)-PD;
PD_error_min = PD-CPSP(1,1);

Y = zeros(4,4);
Y(1,:) = [CR,CR_real,CR_error_max, CR_error_min];
Y(2,:) = [PR,PR_real,PR_error_max, PR_error_min];
Y(3,:) = [SD,SD_real,SD_error_max, SD_error_min];
Y(4,:) = [PD,PD_real,PD_error_max, PD_error_min];
for i =1:4
    for j = 1:4
str= sprintf('%.2g', Y(i,j));
Y(i, j) = str2double(str);
    end
end
name = {'CR','PR','SD','PD'};

figure(1)
for i = 1:4
subplot(1,4,i)
h = bar(Y(i,1:2),'FaceColor', 'flat');
hold on
errorbar(1,Y(i,1),Y(i,3), Y(i,4), 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 2);
hold on
for j = 1:2
text(j-0.25,75,[num2str(Y(i,j)),'%'],'FontSize', 10);
end
set(gca,'YLim',[0,100]);  % Set font name
h.CData(1,:) = cell2mat(color(1,color_size(i)));
h.CData(2,:) = cell2mat(color(7,color_size(i)));
set(gca,'FontSize',14); 
set(gca,'FontName','Airal'); 
ytickformat('percentage');
xticklabels({'S','C'});
title(name(i));
end

set(gcf,'unit','centimeters','position',[1 1 30 3]); 
print('-dpng','-r600','picture\Percentage_change_bar_error.png')

end