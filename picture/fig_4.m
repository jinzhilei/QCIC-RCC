function fig_4()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Plot: This file is used to plot figure 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load data

file='./data/fig.4/CPSP_N.mat';
CPSP_NoDrug=cell2mat(struct2cell(load(file)));
CPSP_NoDrug = fliplr(CPSP_NoDrug);
file='./data/fig.4/CPSP_A.mat';
CPSP_A=cell2mat(struct2cell(load(file)));
CPSP_A = fliplr(CPSP_A);
file='./data/fig.4/CPSP_AB.mat';
CPSP_AB=cell2mat(struct2cell(load(file)));
CPSP_AB = fliplr(CPSP_AB);

file='./data/fig.4/CPSP_N_Ave.mat';
CPSP_NoDrug_ave=cell2mat(struct2cell(load(file)));
CPSP_NoDrug_ave = fliplr(CPSP_NoDrug_ave);

file='./data/fig.4/CPSP_A_Ave.mat';
CPSP_A_ave=cell2mat(struct2cell(load(file)));
CPSP_A_ave = fliplr(CPSP_A_ave);

file='./data/fig.4/CPSP_AB_Ave.mat';
CPSP_AB_ave=cell2mat(struct2cell(load(file)));
CPSP_AB_ave = fliplr(CPSP_AB_ave);

CR(1) = CPSP_NoDrug_ave(1);
PR(1) = CPSP_NoDrug_ave(2);
SD(1) = CPSP_NoDrug_ave(3);
PD(1) = CPSP_NoDrug_ave(4);

CR(2) = CPSP_A_ave(1);
PR(2) = CPSP_A_ave(2);
SD(2) = CPSP_A_ave(3);
PD(2) = CPSP_A_ave(4);

CR(3) = CPSP_AB_ave(1);
PR(3) = CPSP_AB_ave(2);
SD(3) = CPSP_AB_ave(3);
PD(3) = CPSP_AB_ave(4);

%%  real data

CR_real(1) = 0.0;
PR_real(1) = 8.85;
SD_real(1) = 57.1;
PD_real(1) = 34.05;

CR_real(2) = 6.06;
PR_real(2) = 19.19;
SD_real(2) = 36.36;
PD_real(2) = 38.38;

CR_real(3) = 5.66;
PR_real(3) = 33.49;
SD_real(3) = 41.5;
PD_real(3) = 18.47;

CR_error_max(1) = max(CPSP_NoDrug(:,1))-CR(1);
CR_error_min(1) = CR(1)-min(CPSP_NoDrug(:,1));
PR_error_max(1) = max(CPSP_NoDrug(:,2))-PR(1);
PR_error_min(1) = PR(1)-min(CPSP_NoDrug(:,2));
SD_error_max(1) = max(CPSP_NoDrug(:,3))-SD(1);
SD_error_min(1) = SD(1)-min(CPSP_NoDrug(:,3));
PD_error_max(1) = max(CPSP_NoDrug(:,4))-PD(1);
PD_error_min(1) = PD(1)-min(CPSP_NoDrug(:,4));

CR_error_max(2) = max(CPSP_A(:,1))-CR(2);
CR_error_min(2) = CR(2)-min(CPSP_A(:,1));
PR_error_max(2) = max(CPSP_A(:,2))-PR(2);
PR_error_min(2) = PR(2)-min(CPSP_A(:,2));
SD_error_max(2) = max(CPSP_A(:,3))-SD(2);
SD_error_min(2) = SD(2)-min(CPSP_A(:,3));
PD_error_max(2) = max(CPSP_A(:,4))-PD(2);
PD_error_min(2) = PD(2)-min(CPSP_A(:,4));

CR_error_max(3) = max(CPSP_AB(:,1))-CR(3);
CR_error_min(3) = CR(3)-min(CPSP_AB(:,1));
PR_error_max(3) = max(CPSP_AB(:,2))-PR(3);
PR_error_min(3) = PR(3)-min(CPSP_AB(:,2));
SD_error_max(3) = max(CPSP_AB(:,3))-SD(3);
SD_error_min(3) = SD(3)-min(CPSP_AB(:,3));
PD_error_max(3) = max(CPSP_AB(:,4))-PD(3);
PD_error_min(3) = PD(3)-min(CPSP_AB(:,4));

figure(1)

for a = 1:3
    subplot(1,3,a)
    Y = zeros(4,4);
    Y(1,:) = [CR(a),CR_real(a),CR_error_max(a), CR_error_min(a)];
    Y(2,:) = [PR(a),PR_real(a),PR_error_max(a), PR_error_min(a)];
    Y(3,:) = [SD(a),SD_real(a),SD_error_max(a), SD_error_min(a)];
    Y(4,:) = [PD(a),PD_real(a),PD_error_max(a), PD_error_min(a)];
    Y = round(Y,1);
    for i =1:4
        for j = 1:4
            str= sprintf('%.2g', Y(i,j));
            Y(i, j) = str2double(str);
        end
    end
    name = {'CR','PR','SD','PD'};
    T = {'Placebo','Atezolizumab','Atezolizumab+Bevacizumab'};

    h = bar(Y(:,1:2),1,'FaceColor', 'flat');
    hold on
    errorbar([1:4]-0.15,Y(:,1),Y(:,3), Y(:,4), 'LineStyle', 'none', 'Color', [139 137 137]./255, 'LineWidth', 2);
    hold on
    for b =1:4
        text(b-0.28,Y(b,1)+4,[num2str(Y(b,1))],'FontSize', 12);
        text(b+0.03,Y(b,2)+4,[num2str(Y(b,2))],'FontSize', 12);
    end
    set(gca,'YLim',[0,100]);  % Set font name
    h(1).FaceColor = [135 206 250]./255;
    h(2).FaceColor = [240 128 128]./255;
    grid on
    set(gca,'FontSize',12);
    xticklabels(name);
    title(T(a))
    ylabel('Percentage (%)')
    set(gca,'FontName','Airal');
end
subplot(1,3,1)
legend('Model','Data')
legend boxoff
set(gca,'FontName','Airal');

set(gcf,'unit','centimeters','position',[1 1 35 7]);
orient landscape;
print('-dpng','-r200','picture/fig.4.png')
end
