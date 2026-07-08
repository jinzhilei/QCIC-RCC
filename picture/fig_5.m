function fig_5()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Plot: This file is used to plot figure 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lengthe = 35;
hight = 24;
sample = 100;
word_size = 12;
color = ColorMatrix();
day = [1,84:42:360,360];

file='./data/fig.5/ratio_N.mat';
ratio_N=cell2mat(struct2cell(load(file)));
file='./data/fig.5/ratio_A.mat';
ratio_A=cell2mat(struct2cell(load(file)));
file='./data/fig.5/ratio_AB.mat';
ratio_AB=cell2mat(struct2cell(load(file)));
day = [1,84:42:360,360];
ratio_PD = zeros(360,100);
ratio_SD = zeros(360,100);
ratio_PR = zeros(360,100);
ratio_CR = zeros(360,100);
ratio_PD_N = zeros(360,100);
ratio_SD_N = zeros(360,100);
ratio_PR_N = zeros(360,100);
ratio_CR_N = zeros(360,100);
ratio_PD_A = zeros(360,100);
ratio_SD_A = zeros(360,100);
ratio_PR_A = zeros(360,100);
ratio_CR_A = zeros(360,100);
figure(1)
%% No drug
subplot('Position', [0.13, 0.73, 0.21, 0.24]); 
for i = 1:100
     if ratio_N(84,i)>= 20
        ratio_PD_N(:,i) = ratio_N(:,i);
        plot(day,ratio_N(day,i),'-','Color',cell2mat(color(6,1)),'LineWidth',0.2,'DisplayName','PD');
    elseif -30 <= ratio_N(84,i) && ratio_N(84,i) < 20
        ratio_SD_N(:,i) = ratio_N(:,i);
        plot(day,ratio_N(day,i),'-','Color',cell2mat(color(6,4)),'LineWidth',0.5,'DisplayName','SD');
    elseif  -80<= ratio_N(84,i) && ratio_N(84,i) < -30
        ratio_PR_N(:,i) = ratio_N(:,i);
        plot(day,ratio_N(day,i),'-','Color',cell2mat(color(6,10)),'LineWidth',0.5,'DisplayName','PR');
        elseif   ratio_N(84,i) < -80
        ratio_CR_N(:,i) = 0;
        plot(day,ratio_N(day,i),'-','Color',cell2mat(color(6,13)),'LineWidth',0.5,'DisplayName','CR');
     end
hold on
end
ratio_PD_N = ratio_PD_N(:,any(ratio_PD_N,1));
ratio_SD_N = ratio_SD_N(:,any(ratio_SD_N,1));
ratio_PR_N = ratio_PR_N(:,any(ratio_PR_N,1));

ratio_PD_N = mean(ratio_PD_N,2);
ratio_SD_N = mean(ratio_SD_N,2);
ratio_PR_N = mean(ratio_PR_N,2);

plot(day,ratio_PD_N(day),'-','Color',cell2mat(color(3,1)),'LineWidth',3,'DisplayName','PD');
hold on
plot(day,ratio_SD_N(day),'-','Color',cell2mat(color(3,4)),'LineWidth',3,'DisplayName','SD');
hold on
plot(day,ratio_PR_N(day),'-','Color',cell2mat(color(3,10)),'LineWidth',3,'DisplayName','PR');
hold on
plot(day,ratio_PD_N(day),'o','MarkerEdgeColor',[1,1,1],'MarkerFaceColor',cell2mat(color(3,1)),'MarkerSize',8,'LineWidth',1.25,'DisplayName','PD');
hold on
plot(day,ratio_SD_N(day),'o','MarkerEdgeColor',[1,1,1],'MarkerFaceColor',cell2mat(color(3,4)),'MarkerSize',8,'LineWidth',1.25,'DisplayName','SD');
hold on
plot(day,ratio_PR_N(day),'o','MarkerEdgeColor',[1,1,1],'MarkerFaceColor',cell2mat(color(3,10)),'MarkerSize',8,'LineWidth',1.25,'DisplayName','PR');
hold on

grid on

set(gca,'XLim',[0,360]);
set(gca,'YLim',[-100,100]);
set(gca,'Linewidth',1);   % Set outer border
set(gca,'FontSize',word_size);   % Set font size
set(gca,'FontName','Airal');  % Set font name
xlabel('time (day)');
ylabel('Percentage change in tumor (%)'); 
title('Placebo');
%
r_N = ratio_N(84,:);
[~, I_sort] = sort(r_N,'descend');
S_sorted = r_N(I_sort);

subplot('Position', [0.45, 0.73, 0.46, 0.24]); 
b = bar(S_sorted,0.7);
for i =1:100
    b.FaceColor='flat';
    if S_sorted(i)>= 20
        b.CData(i,:) = cell2mat(color(3,1));
    elseif -30 <= S_sorted(i) && S_sorted(i) < 20
        b.CData(i,:) = cell2mat(color(3,4));
    elseif  -80<=S_sorted(i) && S_sorted(i) < -30
        b.CData(i,:) = cell2mat(color(3,10));
        elseif   S_sorted(i) < -80
        b.CData(i,:) = cell2mat(color(3,13));
    end
end
hold on
plot([0,sample+10],[20,20],'--k');
plot([0,sample+10],[-30,-30],'--k');
plot([55.5,55.5],[-100,100],'--k','Color',cell2mat(color(1,11)));
set(gca,'XLim',[0,sample+12]);
set(gca,'YLim',[-100,100]);
set(gca,'Linewidth',1);   % Set outer border
set(gca,'FontSize',word_size);   % Set font size
ylabel('Percentage change in tumor (%)');
title('Placebo');
text(105,60,{'PD'},'FontSize',word_size);
text(105,-5,{'SD'},'FontSize',word_size);
text(101,-50,{'PR/CR'},'FontSize',word_size);
set(gca,'FontName','Airal');  % Set font name
%% A
subplot('Position', [0.13, 0.4, 0.21, 0.24]); 
for i = 1:100
     if ratio_A(84,i)>= 20
        ratio_PD_A(:,i) = ratio_A(:,i);
        plot(day,ratio_A(day,i),'-','Color',cell2mat(color(6,1)),'LineWidth',0.2,'DisplayName','PD');
    elseif -30 <= ratio_A(84,i) && ratio_A(84,i) < 20
        ratio_SD_A(:,i) = ratio_A(:,i);
        plot(day,ratio_A(day,i),'-','Color',cell2mat(color(6,4)),'LineWidth',0.2,'DisplayName','SD');
    elseif  -80<= ratio_A(84,i) && ratio_A(84,i) < -30
        ratio_PR_A(:,i) = ratio_A(:,i);
        plot(day,ratio_A(day,i),'-','Color',cell2mat(color(6,10)),'LineWidth',0.2,'DisplayName','PR');
        elseif   ratio_A(84,i) < -80
        ratio_CR_A(:,i) = ratio_A(:,i);
        plot(day,ratio_A(day,i),'-','Color',cell2mat(color(6,13)),'LineWidth',0.2,'DisplayName','CR');
     end
hold on
end

ratio_PD_A = ratio_PD_A(:,any(ratio_PD_A,1));
ratio_SD_A = ratio_SD_A(:,any(ratio_SD_A,1));
ratio_PR_A = ratio_PR_A(:,any(ratio_PR_A,1));
ratio_CR_A = ratio_CR_A(:,any(ratio_CR_A,1));

ratio_PD_A = mean(ratio_PD_A,2);
ratio_SD_A = mean(ratio_SD_A,2);
ratio_PR_A = mean(ratio_PR_A,2);
ratio_CR_A = mean(ratio_CR_A,2);

plot(day,ratio_PD_A(day),'-','Color',cell2mat(color(3,1)),'LineWidth',3,'DisplayName','PD');
hold on
plot(day,ratio_SD_A(day),'-','Color',cell2mat(color(3,4)),'LineWidth',3,'DisplayName','SD');
hold on
plot(day,ratio_PR_A(day),'-','Color',cell2mat(color(3,10)),'LineWidth',3,'DisplayName','PR');
hold on
plot(day,ratio_PD_A(day),'o','MarkerEdgeColor',[1,1,1],'MarkerFaceColor',cell2mat(color(3,1)),'MarkerSize',8,'LineWidth',1.25,'DisplayName','PD');
hold on
plot(day,ratio_SD_A(day),'o','MarkerEdgeColor',[1,1,1],'MarkerFaceColor',cell2mat(color(3,4)),'MarkerSize',8,'LineWidth',1.25,'DisplayName','SD');
hold on
plot(day,ratio_PR_A(day),'o','MarkerEdgeColor',[1,1,1],'MarkerFaceColor',cell2mat(color(3,10)),'MarkerSize',8,'LineWidth',1.25,'DisplayName','PR');
hold on
grid on

set(gca,'XLim',[0,360]);
set(gca,'YLim',[-100,100]);
set(gca,'Linewidth',1);   % Set outer border
set(gca,'FontSize',word_size);   % Set font size
set(gca,'FontName','Airal');  % Set font name
xlabel('time (day)');
ylabel('Percentage change in tumor (%)'); 
title('Atezolizumab');
%
subplot('Position', [0.45, 0.4, 0.46, 0.24]); 
r_A = ratio_A(84,:);
[~, I_sort] = sort(r_A,'descend');
S_sorted = r_A(I_sort);

b = bar(S_sorted,0.7);
for i =1:100
    b.FaceColor='flat';
    if S_sorted(i)>= 20
        b.CData(i,:) = cell2mat(color(3,1));
    elseif -30 <= S_sorted(i) && S_sorted(i) < 20
        b.CData(i,:) = cell2mat(color(3,4));
    elseif  -80<=S_sorted(i) && S_sorted(i) < -30
        b.CData(i,:) = cell2mat(color(3,10));
        elseif   S_sorted(i) < -80
        b.CData(i,:) = cell2mat(color(3,13));
    end
end
hold on
plot([0,sample+10],[20,20],'--k');
plot([0,sample+10],[-30,-30],'--k');
plot([56.5,56.5],[-100,100],'--k','Color',cell2mat(color(1,11)));
set(gca,'XLim',[0,sample+12]);
set(gca,'YLim',[-100,100]);
set(gca,'Linewidth',1);   % Set outer border
set(gca,'FontSize',word_size);   % Set font size
set(gca,'FontName','Airal');  % Set font name
ylabel('Percentage change in tumor (%)');
title('Atezolizumab');
text(105,60,{'PD'},'FontSize',word_size);
text(105,-5,{'SD'},'FontSize',word_size);
text(101,-50,{'PR/CR'},'FontSize',word_size);
%% A+B
subplot('Position', [0.13, 0.07, 0.21, 0.24]); 
for i = 1:100
     if ratio_AB(84,i)>= 20
        ratio_PD(:,i) = ratio_AB(:,i);
        plot(day,ratio_AB(day,i),'-','Color',cell2mat(color(6,1)),'LineWidth',0.2,'DisplayName','PD');
    elseif -30 <= ratio_AB(84,i) && ratio_AB(84,i) < 20
        ratio_SD(:,i) = ratio_AB(:,i);
        plot(day,ratio_AB(day,i),'-','Color',cell2mat(color(6,4)),'LineWidth',0.2,'DisplayName','SD');
    elseif  -80<= ratio_AB(84,i) && ratio_AB(84,i) < -30
        ratio_PR(:,i) = ratio_AB(:,i);
        plot(day,ratio_AB(day,i),'-','Color',cell2mat(color(6,10)),'LineWidth',0.2,'DisplayName','PR');
        elseif   ratio_AB(84,i) < -80
        ratio_CR(:,i) = ratio_AB(:,i);
        plot(day,ratio_AB(day,i),'-','Color',cell2mat(color(6,13)),'LineWidth',0.2,'DisplayName','CR');
     end
hold on
end

ratio_PD = ratio_PD(:,any(ratio_PD,1));
ratio_SD = ratio_SD(:,any(ratio_SD,1));
ratio_PR = ratio_PR(:,any(ratio_PR,1));
ratio_CR = ratio_CR(:,any(ratio_CR,1));

ratio_PD = mean(ratio_PD,2);
ratio_SD = mean(ratio_SD,2);
ratio_PR = mean(ratio_PR,2);
ratio_CR = mean(ratio_CR,2);

plot(day,ratio_PD(day),'-','Color',cell2mat(color(3,1)),'LineWidth',3,'DisplayName','PD');
hold on
plot(day,ratio_SD(day),'-','Color',cell2mat(color(3,4)),'LineWidth',3,'DisplayName','SD');
hold on
plot(day,ratio_PR(day),'-','Color',cell2mat(color(3,10)),'LineWidth',3,'DisplayName','PR');
hold on
plot(day,ratio_CR(day),'-','Color',cell2mat(color(3,13)),'LineWidth',3,'DisplayName','CR');
hold on
plot(day,ratio_PD(day),'o','MarkerEdgeColor',[1,1,1],'MarkerFaceColor',cell2mat(color(3,1)),'MarkerSize',8,'LineWidth',1.25,'DisplayName','PD');
hold on
plot(day,ratio_SD(day),'o','MarkerEdgeColor',[1,1,1],'MarkerFaceColor',cell2mat(color(3,4)),'MarkerSize',8,'LineWidth',1.25,'DisplayName','SD');
hold on
plot(day,ratio_PR(day),'o','MarkerEdgeColor',[1,1,1],'MarkerFaceColor',cell2mat(color(3,10)),'MarkerSize',8,'LineWidth',1.25,'DisplayName','PR');
hold on
plot(day,ratio_CR(day),'o','MarkerEdgeColor',[1,1,1],'MarkerFaceColor',cell2mat(color(3,13)),'MarkerSize',8,'LineWidth',1.25,'DisplayName','CR');
grid on

set(gca,'XLim',[0,360]);
set(gca,'YLim',[-100,100]);
set(gca,'Linewidth',1);   % Set outer border
set(gca,'FontSize',word_size);   % Set font size
xlabel('time (day)');
ylabel('Percentage change in tumor (%)'); 
title('Atezolizumab + Bevacizumab');
set(gca,'FontName','Airal');  % Set font name
%
subplot('Position', [0.45, 0.07, 0.46, 0.24]); 
ratio_AB = ratio_AB(84,:);
[~, I_sort] = sort(ratio_AB,'descend');
S_sorted = ratio_AB(I_sort);
b = bar(S_sorted,0.7);
for i =1:100
    b.FaceColor='flat';
    if S_sorted(i)>= 20
        b.CData(i,:) = cell2mat(color(3,1));
    elseif -30 <= S_sorted(i) && S_sorted(i) < 20
        b.CData(i,:) = cell2mat(color(3,4));
    elseif  -80<=S_sorted(i) && S_sorted(i) < -30
        b.CData(i,:) = cell2mat(color(3,10));
        elseif   S_sorted(i) < -80
        b.CData(i,:) = cell2mat(color(3,13));
    end
end
hold on
plot([0,sample+10],[20,20],'--k');
plot([0,sample+10],[-30,-30],'--k');
plot([39.5,39.5],[-100,100],'--k','Color',cell2mat(color(1,11)));
set(gca,'XLim',[0,sample+12]);
set(gca,'YLim',[-100,100]);
set(gca,'Linewidth',1);   % Set outer border
set(gca,'FontSize',word_size);   % Set font size
set(gca,'FontName','Airal');  % Set font name
ylabel('Percentage change in tumor (%)');
title('Atezolizumab + Bevacizumab');
text(103,60,{'PD'},'FontSize',word_size);
text(103,-5,{'SD'},'FontSize',word_size);
text(101,-50,{'PR/CR'},'FontSize',word_size);

set(gcf,'unit','centimeters','position',[1 1 lengthe hight]);
print('-dpng','-r600','picture/fig.5.png')

end