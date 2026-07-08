function fig_2()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Plot: This file is used to plot figure 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set parameters

bin_number = 20;
lengthe = 21;
hight = 5;
word_size = 6.5;
fitting_size = 8;
color = ColorMatrix();

file = './data/fig.2/score_JS.mat';
JS=cell2mat(struct2cell(load(file)));
file = ['./data/fig.2/','score_KL.mat'];
KL=cell2mat(struct2cell(load(file)));

file = ['./data/fig.2/','Number_VP_patients.mat'];
number=cell2mat(struct2cell(load(file)));
c= length(number);
B = zeros(c,3);
file=['./data/fig.2/','VP_sample_immune.mat'];
A=cell2mat(struct2cell(load(file)));

CD4 = A(number,12);
CD8 = A(number,21);
Treg = A(number,15);
TAM = A(number,24);
B(:,1) = CD8./CD4;   B(:,2) = CD4./Treg;
B(:,3) = Treg./TAM;

B = log(B);

file = './data/fig.2/True_data_ksdensity.mat';
PDF=cell2mat(struct2cell(load(file)));

vector_data_CD8CD4 = PDF(1,:);
PDF_CD8CD4 = PDF(2,:);
vector_data_CD4Treg = PDF(3,:);
PDF_CD4Treg = PDF(4,:);
vector_data_TregTAM = PDF(5,:);
PDF_TregTAM = PDF(6,:);


%% Plot

figure(1)
subplot(1,3,1)
Histogram_CD8CD4 = histogram(B(:,1),'FaceAlpha',0.5,'Normalization','pdf');
Histogram_CD8CD4.NumBins = bin_number;
hold on
text(-5.25,0.35,['JS: ',num2str(JS(1),'%.4f')],'FontSize',fitting_size,'FontName','Airal','FontWeight','bold');
text(-5.25,0.25,['KL: ',num2str(KL(1),'%.4f')],'FontSize',fitting_size,'FontName','Airal','FontWeight','bold');
Histogram_CD8CD4.FaceColor = cell2mat(color(7,1));
hold on
plot(vector_data_CD8CD4,PDF_CD8CD4,'LineWidth',2,'Color',cell2mat(color(3,1)))
grid on
set(gca,'XLim',[-6 4]);
set(gca,'YLim',[0 0.7]);
set(gca,'Linewidth',1);
set(gca,'FontSize',word_size);
set(gca,'FontName','Airal');
xlabel('Log-transformed CD8/CD4 Ratio');
ylabel('Density');
title('CD8/CD4','FontSize',word_size);
legend('Virtual Patient','Immunogenomic','Location','northwest','FontSize',word_size,'box','off');

subplot(1,3,2)

Histogram_CD4Treg = histogram(B(:,2),'FaceAlpha',0.5,'Normalization','pdf');
Histogram_CD4Treg.NumBins = bin_number;
hold on
text(-4.25,0.35,['JS: ',num2str(JS(2),'%.4f')],'FontSize',fitting_size,'FontName','Airal','FontWeight','bold');
text(-4.25,0.25,['KL: ',num2str(KL(2),'%.4f')],'FontSize',fitting_size,'FontName','Airal','FontWeight','bold');
Histogram_CD4Treg.FaceColor = cell2mat(color(7,1));
hold on
plot(vector_data_CD4Treg,PDF_CD4Treg,'LineWidth',2,'Color',cell2mat(color(3,1)))
grid on
set(gca,'XLim',[-5 4.5]);
set(gca,'YLim',[0 0.7]);
set(gca,'Linewidth',1);
set(gca,'FontSize',word_size);
set(gca,'FontName','Airal');
xlabel('Log-transformed CD4/Treg Ratio');
ylabel('Density');
title('CD4/Treg','FontSize',word_size);

subplot(1,3,3)

Histogram_CD4Treg = histogram(B(:,3),'FaceAlpha',0.5,'Normalization','pdf');
Histogram_CD4Treg.NumBins = bin_number;
hold on
text(-5.25,0.35/0.7*0.6,['JS: ',num2str(JS(3),'%.4f')],'FontSize',fitting_size,'FontName','Airal','FontWeight','bold');
text(-5.25,0.25/0.7*0.6,['KL: ',num2str(KL(3),'%.4f')],'FontSize',fitting_size,'FontName','Airal','FontWeight','bold');
Histogram_CD4Treg.FaceColor = cell2mat(color(7,1));
hold on
plot(vector_data_TregTAM,PDF_TregTAM,'LineWidth',2,'Color',cell2mat(color(3,1)))
grid on
set(gca,'XLim',[-6 4.5]); 
set(gca,'YLim',[0 0.6]);
set(gca,'Linewidth',1);
set(gca,'FontSize',word_size);
set(gca,'FontName','Airal'); 
xlabel('Log-transformed Treg/TAM Ratio');
ylabel('Density');
title('Treg/TAM','FontSize',word_size);

set(gcf,'unit','centimeters','position',[1 1 lengthe hight]);
print('-dpng','-r600','picture\fig.2.png')


end