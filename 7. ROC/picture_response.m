function picture_response()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  plot 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sample = 10000;
number = 100;
day = 43;
lengthe = 15;
hight = 10;
size = 1;
size_tot = 2;
word_size = 18;
cell_sort_AB = zeros(sample,32);
file_response='./data/immune_response.mat';
response = cell2mat(struct2cell(load(file_response)));
color = ColorMatrix();
color_size = [1,2,3,4,7,8,10,11,12,14];

for i =1:sample
    file=['../VP_AB/VP_AB_bar_number/','VPsample_AB',num2str(i),'.mat'];
    cell=cell2mat(struct2cell(load(file)));
    cell_sort_AB(i,:) = cell(day,:);
end
%Th
for j = 1:20
L = randi([1,sample],number,1);
[A, B, ~, ~] = perfcurve(response(3,L), cell_sort_AB(L,12), 0);  
actual_points = length(A);
    

    X(1:actual_points,j) = A;
    Y(1:actual_points,j) = B;
end

[C(:,1), D(:,1), ~, AUC_t] = perfcurve(response(3,:), cell_sort_AB(:,12), 0);  
figure(1)
for i = 1:20
str= sprintf('%.2g', AUC_t);
AUC_t = str2double(str);
plot(X(:,i), Y(:,i), '-', 'LineWidth', size,'Color',cell2mat(color(7,color_size(1))));  % X=FPR, Y=TPR
hold on;
end
plot(C(:,1), D(:,1), '-', 'LineWidth',size_tot ,'Color',cell2mat(color(1,color_size(1))));  % X=FPR, Y=TPR
hold on
plot([0 1], [0 1], '--','Color',[0.2,0.2,0.2],'LineWidth', size_tot);  
% text(0.8,0.2,'AUC =' num2str(AUC_t),'FontSize',1);
xlabel('1-Specificity');
ylabel('Sensitivity');
title('Atezolizumab + Bevacizumab')
grid on;
legend(['Th           AUC =' num2str(AUC_t)], 'Location', 'southeast');
set(gca,'FontSize',word_size);   % Set font size
set(gcf,'unit','centimeters','position',[1 1 lengthe hight]); 
print('-dpng','-r600','picture\responce_AB_Th.png')
% Tr
for j = 1:20
L = randi([1,sample],number,1);
[A, B, ~, ~] = perfcurve(response(3,L), cell_sort_AB(L,15), 1);  
actual_points = length(A);
    
    X(1:actual_points,j) = A;
    Y(1:actual_points,j) = B;
end

[C, D, ~, AUC_t] = perfcurve(response(3,:), cell_sort_AB(:,15), 1); 

figure(2)
for i = 1:20
str= sprintf('%.2g', AUC_t);
AUC_t = str2double(str);
plot(X(:,i), Y(:,i), '-', 'LineWidth', size,'Color',cell2mat(color(7,color_size(2))));  % X=FPR, Y=TPR
hold on;
end
plot(C, D, '-', 'LineWidth',size_tot ,'Color',cell2mat(color(1,color_size(2))));  % X=FPR, Y=TPR
hold on
plot([0 1], [0 1], '--','Color',[0.2,0.2,0.2],'LineWidth', size_tot);  
% text(0.8,0.2,'AUC =' num2str(AUC_t),'FontSize',1);
xlabel('1-Specificity');
ylabel('Sensitivity');
title('Atezolizumab + Bevacizumab')
grid on;
legend(['Tr           AUC =' num2str(AUC_t)], 'Location', 'southeast');
set(gca,'FontSize',word_size);   % Set font size
set(gcf,'unit','centimeters','position',[1 1 lengthe hight]); 
print('-dpng','-r600','picture\responce_AB_Tr.png')

% Tc
for j = 1:20
L = randi([1,sample],number,1);
[A, B, ~, ~] = perfcurve(response(3,L), cell_sort_AB(L,21), 1);  
actual_points = length(A);
    

    X(1:actual_points,j) = A;
    Y(1:actual_points,j) = B;
end

[C, D, ~, AUC_t] = perfcurve(response(3,:), cell_sort_AB(:,21), 1); 

figure(3)
for i = 1:20
str= sprintf('%.2g', AUC_t);
AUC_t = str2double(str);
plot(X(:,i), Y(:,i), '-', 'LineWidth', size,'Color',cell2mat(color(7,color_size(3))));  % X=FPR, Y=TPR
hold on;
end
plot(C, D, '-', 'LineWidth',size_tot ,'Color',cell2mat(color(1,color_size(4))));  % X=FPR, Y=TPR
hold on
plot([0 1], [0 1], '--','Color',[0.2,0.2,0.2],'LineWidth', size_tot);  
% text(0.8,0.2,'AUC =' num2str(AUC_t),'FontSize',1);
xlabel('1-Specificity');
ylabel('Sensitivity');
title('Atezolizumab + Bevacizumab')
grid on;
legend(['Tc           AUC =' num2str(AUC_t)], 'Location', 'southeast');
set(gca,'FontSize',word_size);   % Set font size
set(gca,'FontSize',word_size);   % Set font size
set(gcf,'unit','centimeters','position',[1 1 lengthe hight]); 
print('-dpng','-r600','picture\responce_AB_Tc.png')

% TAM
for j = 1:20
L = randi([1,sample],number,1);
[A, B, ~, ~] = perfcurve(response(3,L), cell_sort_AB(L,24), 1);  
actual_points = length(A);
    
    X(1:actual_points,j) = A;
    Y(1:actual_points,j) = B;
end

[C, D, ~, AUC_t] = perfcurve(response(3,:), cell_sort_AB(:,24), 1);  

figure(4)
for i = 1:20
str= sprintf('%.2g', AUC_t);
AUC_t = str2double(str);
plot(X(:,i), Y(:,i), '-', 'LineWidth', size,'Color',cell2mat(color(7,color_size(4))));  % X=FPR, Y=TPR
hold on;
end
plot(C, D, '-', 'LineWidth',size_tot ,'Color',cell2mat(color(1,color_size(4))));  % X=FPR, Y=TPR
hold on
plot([0 1], [0 1], '--','Color',[0.2,0.2,0.2],'LineWidth', size_tot); 
% text(0.8,0.2,'AUC =' num2str(AUC_t),'FontSize',1);
xlabel('1-Specificity');
ylabel('Sensitivity');
title('Atezolizumab + Bevacizumab')
grid on;
legend(['TAM           AUC =' num2str(AUC_t)], 'Location', 'southeast');
set(gca,'FontSize',word_size);   % Set font size
set(gcf,'unit','centimeters','position',[1 1 lengthe hight]); 
print('-dpng','-r600','picture\responce_AB_TAM.png')

% Th/Tr
for j = 1:20
L = randi([1,sample],number,1);
[A, B, ~, ~] = perfcurve(response(3,L), cell_sort_AB(L,12)./cell_sort_AB(L,15), 0);  
actual_points = length(A);
    X(1:actual_points,j) = A;
    Y(1:actual_points,j) = B;
end
[C, D, ~, AUC_t] = perfcurve(response(3,:), cell_sort_AB(:,12)./cell_sort_AB(:,15), 0);  
figure(5)
for i = 1:20
str= sprintf('%.2g', AUC_t);
AUC_t = str2double(str);
plot(X(:,i), Y(:,i), '-', 'LineWidth', size,'Color',cell2mat(color(7,color_size(5))));  % X=FPR, Y=TPR
hold on;
end
plot(C, D, '-', 'LineWidth',size_tot ,'Color',cell2mat(color(1,color_size(5))));  % X=FPR, Y=TPR
hold on
plot([0 1], [0 1], '--','Color',[0.2,0.2,0.2],'LineWidth', size_tot);  
% text(0.8,0.2,'AUC =' num2str(AUC_t),'FontSize',1);
xlabel('1-Specificity');
ylabel('Sensitivity');
title('Atezolizumab + Bevacizumab')
grid on;
legend(['Th/Tr           AUC =' num2str(AUC_t)], 'Location', 'southeast');
set(gca,'FontSize',word_size);   % Set font size
set(gcf,'unit','centimeters','position',[1 1 lengthe hight]); 
print('-dpng','-r600','picture\responce_AB_Th_Tr.png')

% Th/Tc
for j = 1:20
L = randi([1,sample],number,1);
[A, B, ~, ~] = perfcurve(response(3,L), cell_sort_AB(L,12)./cell_sort_AB(L,21), 0);  
actual_points = length(A);
    X(1:actual_points,j) = A;
    Y(1:actual_points,j) = B;
end
[C, D, ~, AUC_t] = perfcurve(response(3,:), cell_sort_AB(:,12)./cell_sort_AB(:,21), 0);  
figure(6)
for i = 1:20
str= sprintf('%.2g', AUC_t);
AUC_t = str2double(str);
plot(X(:,i), Y(:,i), '-', 'LineWidth', size,'Color',cell2mat(color(7,color_size(6))));  % X=FPR, Y=TPR
hold on;
end
plot(C, D, '-', 'LineWidth',size_tot ,'Color',cell2mat(color(1,color_size(6))));  % X=FPR, Y=TPR
hold on
plot([0 1], [0 1], '--','Color',[0.2,0.2,0.2],'LineWidth', size_tot);  
% text(0.8,0.2,'AUC =' num2str(AUC_t),'FontSize',1);
xlabel('1-Specificity');
ylabel('Sensitivity');
title('Atezolizumab + Bevacizumab')
grid on;
legend(['Th/Tc           AUC =' num2str(AUC_t)], 'Location', 'southeast');
set(gca,'FontSize',word_size);   % Set font size
set(gcf,'unit','centimeters','position',[1 1 lengthe hight]); 
print('-dpng','-r600','picture\responce_AB_Th_Tc.png')

% Th/TAM
for j = 1:20
L = randi([1,sample],number,1);
[A, B, ~, ~] = perfcurve(response(3,L), cell_sort_AB(L,12)./cell_sort_AB(L,24), 0);  
actual_points = length(A);
    X(1:actual_points,j) = A;
    Y(1:actual_points,j) = B;
end
[C, D, ~, AUC_t] = perfcurve(response(3,:), cell_sort_AB(:,12)./cell_sort_AB(:,24), 0); 
figure(7)
for i = 1:20
str= sprintf('%.2g', AUC_t);
AUC_t = str2double(str);
plot(X(:,i), Y(:,i), '-', 'LineWidth', size,'Color',cell2mat(color(7,color_size(7))));  % X=FPR, Y=TPR
hold on;
end
plot(C, D, '-', 'LineWidth', size_tot ,'Color',cell2mat(color(1,color_size(7))));  % X=FPR, Y=TPR
hold on
plot([0 1], [0 1], '--','Color',[0.2,0.2,0.2],'LineWidth', size_tot);  
% text(0.8,0.2,'AUC =' num2str(AUC_t),'FontSize',1);
xlabel('1-Specificity');
ylabel('Sensitivity');
title('Atezolizumab + Bevacizumab')
grid on;
legend(['Th/TAM           AUC =' num2str(AUC_t)], 'Location', 'southeast');
set(gca,'FontSize',word_size);   % Set font size
set(gcf,'unit','centimeters','position',[1 1 lengthe hight]); 
print('-dpng','-r600','picture\responce_AB_Th_TAM.png')

% Tr/Tc
for j = 1:20
L = randi([1,sample],number,1);
[A, B, ~, ~] = perfcurve(response(3,L), cell_sort_AB(L,15)./cell_sort_AB(L,21), 0);  
actual_points = length(A);
    X(1:actual_points,j) = A;
    Y(1:actual_points,j) = B;
end
[C, D, ~, AUC_t] = perfcurve(response(3,:), cell_sort_AB(:,15)./cell_sort_AB(:,21), 0);  
figure(8)
for i = 1:20
str= sprintf('%.2g', AUC_t);
AUC_t = str2double(str);
plot(X(:,i), Y(:,i), '-', 'LineWidth', size,'Color',cell2mat(color(7,color_size(8))));  % X=FPR, Y=TPR
hold on;
end
plot(C, D, '-', 'LineWidth', size_tot ,'Color',cell2mat(color(1,color_size(8))));  % X=FPR, Y=TPR
hold on
plot([0 1], [0 1], '--','Color',[0.2,0.2,0.2],'LineWidth', size_tot);  
% text(0.8,0.2,'AUC =' num2str(AUC_t),'FontSize',1);
xlabel('1-Specificity');
ylabel('Sensitivity');
title('Atezolizumab + Bevacizumab')
grid on;
legend(['Tr/Tc           AUC =' num2str(AUC_t)], 'Location', 'southeast');
set(gca,'FontSize',word_size);   % Set font size
set(gcf,'unit','centimeters','position',[1 1 lengthe hight]); 
print('-dpng','-r600','picture\responce_AB_Tr_Tc.png')

% Tr/TAM
for j = 1:20
L = randi([1,sample],number,1);
[A, B, ~, ~] = perfcurve(response(3,L), cell_sort_AB(L,15)./cell_sort_AB(L,24), 1);  
actual_points = length(A);

    X(1:actual_points,j) = A;
    Y(1:actual_points,j) = B;
end
[C, D, ~, AUC_t] = perfcurve(response(3,:), cell_sort_AB(:,15)./cell_sort_AB(:,24), 1);  
figure(9)
for i = 1:20
str= sprintf('%.2g', AUC_t);
AUC_t = str2double(str);
plot(X(:,i), Y(:,i), '-', 'LineWidth', size,'Color',cell2mat(color(7,color_size(9))));  % X=FPR, Y=TPR
hold on;
end
plot(C, D, '-', 'LineWidth', size_tot ,'Color',cell2mat(color(1,color_size(9))));  % X=FPR, Y=TPR
hold on
plot([0 1], [0 1], '--','Color',[0.2,0.2,0.2],'LineWidth', size_tot);  
% text(0.8,0.2,'AUC =' num2str(AUC_t),'FontSize',1);
xlabel('1-Specificity');
ylabel('Sensitivity');
title('Atezolizumab + Bevacizumab')
grid on;
legend(['Tr/TAM           AUC =' num2str(AUC_t)], 'Location', 'southeast');
set(gca,'FontSize',word_size);   % Set font size
set(gcf,'unit','centimeters','position',[1 1 lengthe hight]); 
print('-dpng','-r600','picture\responce_AB_Tr_TAM.png')

% Tc/TAM
for j = 1:20
L = randi([1,sample],number,1);
[A, B, ~, ~] = perfcurve(response(3,L), cell_sort_AB(L,21)./cell_sort_AB(L,24), 1);  
actual_points = length(A);
    X(1:actual_points,j) = A;
    Y(1:actual_points,j) = B;
end
[C, D, ~, AUC_t] = perfcurve(response(3,:), cell_sort_AB(:,21)./cell_sort_AB(:,24), 1);  
figure(10)
for i = 1:20
str= sprintf('%.2g', AUC_t);
AUC_t = str2double(str);
plot(X(:,i), Y(:,i), '-', 'LineWidth', size,'Color',cell2mat(color(7,color_size(10))));  % X=FPR, Y=TPR
hold on;
end
plot(C, D, '-', 'LineWidth', size_tot,'Color',cell2mat(color(1,color_size(10))));  % X=FPR, Y=TPR
hold on
plot([0 1], [0 1], '--','Color',[0.2,0.2,0.2],'LineWidth', size_tot);  
xlabel('1-Specificity');
ylabel('Sensitivity');
title('Atezolizumab + Bevacizumab')
grid on;
legend(['Tc/TAM           AUC =' num2str(AUC_t)], 'Location', 'southeast');
set(gca,'FontSize',word_size);   % Set font size
set(gcf,'unit','centimeters','position',[1 1 lengthe hight]); 
print('-dpng','-r600','picture\responce_AB_Tc_TAM.png')

end

