function picture_all_recurrence_time()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% This function is used to plot the recurrence time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

number = 192;
lengthe = 50;
hight = 50;
xlim = 200;
ylim = 220;
word_size = 15;
A = zeros(10000,number);
C = zeros(1,number);
Z = zeros(24,5);
value = -1;
Y = {'(1,1)','(1,2)','(1,3)','(1,4)','(1,5)','(1,6)','(1,7)','(1,8)','(2,1)','(2,2)','(2,3)','(2,4)','(2,5)','(2,6)','(2,7)','(2,8)',...
    '(3,1)','(3,2)','(3,3)','(3,4)','(3,5)','(3,6)','(3,7)','(3,8)','(4,1)','(4,2)','(4,3)','(4,4)','(4,5)','(4,6)','(4,7)','(4,8)',...
    '(5,1)','(5,2)','(5,3)','(5,4)','(5,5)','(5,6)','(5,7)','(5,8)','(6,1)','(6,2)','(6,3)','(6,4)','(6,5)','(6,6)','(6,7)','(6,8)',...
    '(7,1)','(7,2)','(7,3)','(7,4)','(7,5)','(7,6)','(7,7)','(7,8)','(8,1)','(8,2)','(8,3)','(8,4)','(8,5)','(8,6)','(8,7)','(8,8)',};
for p = 1:3
    for m = 1:8
        for n = 1:8
            file=strcat('Recurrence_rate/treat',num2str(p),'_Recurrence_rate(',num2str(m),',',num2str(n),').mat');
            A(:,m+(n-1)*8+(p-1)*64)=cell2mat(struct2cell(load(file)));
        end
    end
end

for i = 1:number
B = A(:,i);
B(B==value) = [];
a = mean(B);
a = floor(a);
C(i) = a;
end

figure(1)
subplot(1,3,1)
scatter(C(1:64),1:64)
hold on
for i =1:64
plot([0,C(i)],[i,i],'-k')
hold on
end
a= mean(C(1:64));
a = floor(a);
plot([a,a],[1-5,69],'--')
hold on
h = text(a+2,30,['$\textrm{mean: }$',num2str(a),'$\textrm{day}$'],"Interpreter","latex");
h.FontSize = word_size;            
set(gca,'XLim',[0 ylim]); 
set(gca,'YLim',[1-5,69]);  % Set coordinate range (Y)
set(gca, 'YTick', 1:64);  
set(gca, 'YTickLabel', Y); 
set(gca, 'FontName', 'Airal');  
title('Clinical administration');
xlabel('Time (day)','FontSize',word_size)


subplot(1,3,2)
scatter(C(65:128),65:128)
hold on
for i =65:128
plot([0,C(i)],[i,i],'-k')
hold on
end
a= mean(C(65:128));
a = floor(a);
plot([a,a],[60,133],'--')
hold on
h = text(a+2,94,['$\textrm{mean: }$',num2str(a),'$\textrm{day}$'],"Interpreter","latex");
h.FontSize = word_size;           
set(gca,'XLim',[0 ylim]);  
set(gca,'YLim',[60,133]);  % Set coordinate range (Y)
set(gca, 'YTick', 65:128); 
set(gca, 'YTickLabel', Y);  
title('Sustained dosing');
xlabel('Time (day)','FontSize',word_size)
set(gca, 'FontName', 'Airal'); 

subplot(1,3,3)
scatter(C(129:number),129:number)
hold on
for i =129:number
plot([0,C(i)],[i,i],'-k')
hold on
end
a= mean(C(129:number));
a = floor(a);
plot([a,a],[129-5,number+5],'--')
hold on
h = text(a+2,160-2,['$\textrm{mean: }$',num2str(a),'$\textrm{day}$'],"Interpreter","latex");
h.FontSize = word_size;              
set(gca,'XLim',[0 ylim]);  
set(gca,'YLim',[129-5,number+5]);  % Set coordinate range (Y)
set(gca, 'YTick', 129:number);  
set(gca, 'YTickLabel', Y);  
title('Waning treatment');
xlabel('Time (day)','FontSize',word_size)
set(gca, 'FontName', 'Airal'); 

set(gcf,'unit','centimeters','position',[1 1 lengthe hight]); 
print('-dpng','-r600','picture/recurrence_all_time.png')

for p =1:3
for m = 1:8
    for n = 1:8
     Z((p-1)*8+m,n) = C(n+(m-1)*8+(p-1)*64);
    end
end
end

figure(2)
m = 1:8;
n = 1:8;
h = surfc(m,n,Z(1:8,:));
hold on
[C,p]=contour(Z(1:8,:)); 
clabel(C,p); 
ylabel('Atezolizumab');
xlabel('Bevacizumab');
title('Clinical administration');
h(1).FaceColor = 'interp';
set(gca,'FontSize',word_size);   % Set font size
set(gca,'FontName','Airal');  % Set font name
set(gca,'ZLim',[0 xlim]);  % Set coordinate range (X)


figure(3)
m = 1:8;
n = 1:8;
h = surfc(m,n,Z(9:16,:));
hold on
[C,p]=contour(Z(9:16,:)); 
clabel(C,p); 
ylabel('Atezolizumab');
xlabel('Bevacizumab');
title('Sustained dosing');
h(1).FaceColor = 'interp';
set(gca,'FontSize',word_size);   % Set font size
set(gca,'FontName','Airal');  % Set font name
set(gca,'ZLim',[0 xlim]);  % Set coordinate range (X)

% 
figure(4)
m = 1:8;
n = 1:8;
h = surfc(m,n,Z(17:24,:));
hold on
[C,p]=contour(Z(17:24,:)); 
clabel(C,p); 
xlabel('Atezolizumab');
ylabel('Bevacizumab');
title('Waning treatment');
h(1).FaceColor = 'interp';
set(gca,'FontSize',word_size);   % Set font size
set(gca,'FontName','Airal');  % Set font name
set(gca,'ZLim',[0 xlim]);  % Set coordinate range (X)


end