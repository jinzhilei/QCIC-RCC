function picture_bar_error()
color = ColorMatrix();
color_size = [7,10,4,1];

file='./VP_NoDrug_bar/CPSP.mat';
CPSP=cell2mat(struct2cell(load(file)));
CPSP = sort(CPSP);
CR = mean(CPSP(:,4));
PR = mean(CPSP(:,3));
SD = mean(CPSP(:,2));
PD = mean(CPSP(:,1));

CR_real = 0;
PR_real = 8.85;
SD_real = 57.1;
PD_real = 34.05;

CR_error_max = CPSP(100,4)-CR;
CR_error_min = CR-CPSP(1,4);
PR_error_max = CPSP(100,3)-PR;
PR_error_min = PR-CPSP(1,3);
SD_error_max = CPSP(100,2)-SD;
SD_error_min = SD-CPSP(1,2);
PD_error_max = CPSP(100,1)-PD;
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
        text(j-0.25,70,[num2str(Y(i,j)),'%'],'FontSize', 14);
    end

    set(gca,'YLim',[0,80]);  % Set font name
    h.CData(1,:) = cell2mat(color(1,color_size(i)));
    h.CData(2,:) = cell2mat(color(7,color_size(i)));
    set(gca,'FontSize',14);
    set(gca,'FontName','Airal');
    ytickformat('percentage');
    xticklabels({'S','C'});
    title(name(i));

end

set(gcf,'unit','centimeters','position',[1 1 45 5]);
print('-dpng','-r600','picture\Percentage_change_bar_reeor.png')

end