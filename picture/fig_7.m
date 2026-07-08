function fig_7()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  plot 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sample = 1664;
number = 100;
lengthe =40;
hight = 15;
size = 1;
size_tot = 3;
word_size = 13;
point_size = 13;
AUC = zeros(1,7);
AUC_T = zeros(10,7);
rng(42)
file_response='./data/fig.7/response.mat';
response = cell2mat(struct2cell(load(file_response)));
file_response='./data/fig.7/immune.mat';
cell_sort_AB = cell2mat(struct2cell(load(file_response)));
Xlabel = {'Th','Tr','Tc','TAM','Th/Tc','Tr/Tc','Tc/TAM'};
color_box = [1.00, 0.70, 0.70;
             1.00, 0.80, 0.55;
             0.98, 0.92, 0.65;
             0.70, 0.90, 0.70;
             0.65, 0.85, 0.90;
             0.75, 0.70, 0.95;
             0.95, 0.70, 0.85];
softness = 0.4;  
original = [1.00, 0.70, 0.70;
             1.00, 0.80, 0.55;
             0.98, 0.92, 0.65;
             0.70, 0.90, 0.70;
             0.65, 0.85, 0.90;
             0.75, 0.70, 0.95;
             0.95, 0.70, 0.85];
color_piont = original * (1 - softness) + softness * 1;
color_p = [1.00, 0.70, 0.70;
             1.00, 0.80, 0.55;
             0.98, 0.92, 0.65;
             0.70, 0.90, 0.70;
             0.65, 0.85, 0.90;
             0.75, 0.70, 0.95;
             0.95, 0.70, 0.85] *0.95;

%Th
for j = 1:10
L = randi([1,sample],number,1);
[A, B, ~, U] = perfcurve(response(L), cell_sort_AB(L,12), 0);  
actual_points = length(A);
    X(1:actual_points,j) = A;
    Y(1:actual_points,j) = B;
    AUC_T(j,1) = U;
end
[C(:,1), D(:,1), ~, AUC_t] = perfcurve(response, cell_sort_AB(:,12), 0);  
AUC(1) = AUC_t;
figure(1)
subplot(2,5,7)
for i = 1:10
str= sprintf('%.2g', AUC_t);
AUC_t = str2double(str);
plot(X(:,i), Y(:,i), '-', 'LineWidth', size,'Color',color_box(7,:));  % X=FPR, Y=TPR
hold on;
end
plot(C(:,1), D(:,1), '-', 'LineWidth',size_tot ,'Color',color_p(7,:));  % X=FPR, Y=TPR
hold on
plot([0 1], [0 1], '--','Color',[0.8,0.8,0.8],'LineWidth', size_tot-1.5);  
xlabel('1-Specificity');
ylabel('Sensitivity');
grid on;
title('Th ','FontSize',word_size);
set(gca,'FontSize',word_size);   % Set font size

% Tr
for j = 1:10
L = randi([1,sample],number,1);
[A, B, ~, U] = perfcurve(response(L), cell_sort_AB(L,15), 0);  
actual_points = length(A);
    X(1:actual_points,j) = A;
    Y(1:actual_points,j) = B;
    AUC_T(j,2) = U;
end
[C, D, ~, AUC_t] = perfcurve(response, cell_sort_AB(:,15), 0); 
AUC(2) = AUC_t;
subplot(2,5,5)
for i = 1:10
str= sprintf('%.2g', AUC_t);
AUC_t = str2double(str);
plot(X(:,i), Y(:,i), '-', 'LineWidth', size,'Color',color_box(5,:));  % X=FPR, Y=TPR
hold on;
end
plot(C, D, '-', 'LineWidth',size_tot ,'Color',color_p(5,:));  % X=FPR, Y=TPR
hold on
plot([0 1], [0 1], '--','Color',[0.8,0.8,0.8],'LineWidth', size_tot-1.5);  
xlabel('1-Specificity');
ylabel('Sensitivity');
grid on;
title('Tr ');
set(gca,'FontSize',word_size);   % Set font size

% Tc
for j = 1:10
L = randi([1,sample],number,1);
[A, B, ~, U] = perfcurve(response(L), cell_sort_AB(L,21), 1);  
actual_points = length(A);
    X(1:actual_points,j) = A;
    Y(1:actual_points,j) = B;
    AUC_T(j,3) = U;
end
[C, D, ~, AUC_t] = perfcurve(response, cell_sort_AB(:,21), 1); 
AUC(3) = AUC_t;
subplot(2,5,1)
for i = 1:10
str= sprintf('%.2g', AUC_t);
AUC_t = str2double(str);
plot(X(:,i), Y(:,i), '-', 'LineWidth', size,'Color',color_box(1,:));  % X=FPR, Y=TPR
hold on;
end
plot(C, D, '-', 'LineWidth',size_tot ,'Color',color_p(1,:));  % X=FPR, Y=TPR
hold on
plot([0 1], [0 1], '--','Color',[0.8,0.8,0.8],'LineWidth', size_tot-1.5);  
xlabel('1-Specificity');
ylabel('Sensitivity');
grid on;
title('Tc ');
set(gca,'FontSize',word_size);   % Set font size

% TAM
for j = 1:10
L = randi([1,sample],number,1);
[A, B, ~, U] = perfcurve(response(L), cell_sort_AB(L,24), 1);  
actual_points = length(A);
    X(1:actual_points,j) = A;
    Y(1:actual_points,j) = B;
    AUC_T(j,4) = U;
end
[C, D, ~, AUC_t] = perfcurve(response, cell_sort_AB(:,24), 1);  
AUC(4) = AUC_t;
subplot(2,5,6)
for i = 1:10
str= sprintf('%.2g', AUC_t);
AUC_t = str2double(str);
plot(X(:,i), Y(:,i), '-', 'LineWidth', size,'Color',color_box(6,:));  % X=FPR, Y=TPR
hold on;
end
plot(C, D, '-', 'LineWidth',size_tot ,'Color',color_p(6,:));  % X=FPR, Y=TPR
hold on
plot([0 1], [0 1], '--','Color',[0.8,0.8,0.8],'LineWidth', size_tot-1.5); 
xlabel('1-Specificity');
ylabel('Sensitivity');
grid on;
title('TAM ');
set(gca,'FontSize',word_size);   % Set font size

% Th/Tc
for j = 1:10
L = randi([1,sample],number,1);
[A, B, ~, U] = perfcurve(response(L), cell_sort_AB(L,12)./cell_sort_AB(L,21), 0);  
actual_points = length(A);
    X(1:actual_points,j) = A;
    Y(1:actual_points,j) = B;
    AUC_T(j,5) = U;
end
[C, D, ~, AUC_t] = perfcurve(response, cell_sort_AB(:,12)./cell_sort_AB(:,21), 0);  
AUC(5) = AUC_t;
subplot(2,5,4)
for i = 1:10
str= sprintf('%.2g', AUC_t);
AUC_t = str2double(str);
plot(X(:,i), Y(:,i), '-', 'LineWidth', size,'Color',color_box(4,:));  % X=FPR, Y=TPR
hold on;
end
plot(C, D, '-', 'LineWidth',size_tot ,'Color',color_p(4,:));  % X=FPR, Y=TPR
hold on
plot([0 1], [0 1], '--','Color',[0.8,0.8,0.8],'LineWidth', size_tot-1.5);  
xlabel('1-Specificity');
ylabel('Sensitivity');
 
grid on;
title('Th/Tc ');
set(gca,'FontSize',word_size);   % Set font size

% Tr/Tc
for j = 1:10
L = randi([1,sample],number,1);
[A, B, ~, U] = perfcurve(response(L), cell_sort_AB(L,15)./cell_sort_AB(L,21), 0);  
actual_points = length(A);
    X(1:actual_points,j) = A;
    Y(1:actual_points,j) = B;
    AUC_T(j,6) = U;
end
[C, D, ~, AUC_t] = perfcurve(response, cell_sort_AB(:,15)./cell_sort_AB(:,21), 0);  
AUC(6) = AUC_t;
subplot(2,5,2)
for i = 1:10
str= sprintf('%.2g', AUC_t);
AUC_t = str2double(str);
plot(X(:,i), Y(:,i), '-', 'LineWidth', size,'Color',color_box(2,:));  % X=FPR, Y=TPR
hold on;
end
plot(C, D, '-', 'LineWidth', size_tot ,'Color',color_p(2,:));  % X=FPR, Y=TPR
hold on
plot([0 1], [0 1], '--','Color',[0.8,0.8,0.8],'LineWidth', size_tot-1.5);  
xlabel('1-Specificity');
ylabel('Sensitivity');
grid on;
title('Tr/Tc ');
set(gca,'FontSize',word_size);   % Set font size

% Tc/TAM
for j = 1:10
L = randi([1,sample],number,1);
[A, B, ~, U] = perfcurve(response(L), cell_sort_AB(L,21)./cell_sort_AB(L,24), 1);  
actual_points = length(A);
    X(1:actual_points,j) = A;
    Y(1:actual_points,j) = B;
    AUC_T(j,7) = U;
end
[C, D, ~, AUC_t] = perfcurve(response, cell_sort_AB(:,21)./cell_sort_AB(:,24), 1); 
AUC(7) = AUC_t;
subplot(2,5,3)
for i = 1:10
str= sprintf('%.2g', AUC_t);
AUC_t = str2double(str);
plot(X(:,i), Y(:,i), '-', 'LineWidth', size,'Color',color_box(3,:));  % X=FPR, Y=TPR
hold on;
end
plot(C, D, '-', 'LineWidth', size_tot,'Color',color_p(3,:));  % X=FPR, Y=TPR
hold on
plot([0 1], [0 1], '--','Color',[0.8,0.8,0.8],'LineWidth', size_tot-1.5);  
xlabel('1-Specificity');
ylabel('Sensitivity');
grid on;
title('Tc/TAM ');
set(gca,'FontSize',word_size);   % Set font size


subplot('Position', [0.48, 0.12, 0.43, 0.32]); 
[~, I_sort] = sort(AUC,'descend');
S_sorted = AUC_T(:,I_sort); 
x=0.775+rand(1,10)*0.45;
boxplot(S_sorted);
hBox = findobj(gca, 'Tag', 'Box');
c_P = flipud(color_piont);
for i = 1:length(hBox)
    set(hBox(i), 'Color', c_P(i,:), 'LineWidth', 1.5);
    xd = get(hBox(i), 'XData');
    yd = get(hBox(i), 'YData');
    patch(xd, yd, c_P(i,:), 'FaceAlpha', 0.8, 'EdgeColor', 'none');
end
hold on
meanVals = median(AUC(I_sort), 1); 
xPos = 1:7;
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
plot(x,S_sorted(:,1),'.','Color',color_p(1,:),'MarkerSize', point_size);
hold on
plot(x+1,S_sorted(:,2),'.','Color',color_p(2,:),'MarkerSize', point_size)
hold on
plot(x+2,S_sorted(:,3),'.','Color',color_p(3,:),'MarkerSize', point_size)
hold on
plot(x+3,S_sorted(:,4),'.','Color',color_p(4,:),'MarkerSize', point_size)
hold on
plot(x+4,S_sorted(:,5),'.','Color',color_p(5,:),'MarkerSize', point_size)
hold on
plot(x+5,S_sorted(:,6),'.','Color',color_p(6,:),'MarkerSize', point_size)
hold on
plot(x+6,S_sorted(:,7),'.','Color',color_p(7,:),'MarkerSize', point_size)
grid on

xticks(1:7);
ax = gca;
ax.XTickLabel=Xlabel;
ax.XTickLabel = ax.XTickLabel(I_sort);
ylabel('AUC');
ax.LineWidth=1;
ax.FontSize=word_size;
ax.FontName='Airal';
ylim([0.4,0.9])

set(gcf,'unit','centimeters','position',[1 1 lengthe hight]); 
print('-dpng','-r600','picture\fig.7.png')



end

