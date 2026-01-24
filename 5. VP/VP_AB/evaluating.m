function evaluating()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This module is used to describe the tumor evaluation metrics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter sampling

day = 360;
jianc =  43;
T = 21;
lengthe = 8;
hight = 6;
word_size = 9;
L = 0.5;
sample = 100;
tumor = zeros(day,sample);
ratio = zeros(day,sample);
A = zeros(1,4);
color = ColorMatrix();

file_R='./VP_AB/VP_tumor_number.mat';
r=cell2mat(struct2cell(load(file_R)));
tumor_0 = ((r.*2.572e-9.*0.8*3)./(0.37*4*3.14)).^(1/3);

%% Calculation of tumor evaluation indicators

for i =1:sample
    file=['./VP_AB/','VPsample_AB',num2str(i),'.mat'];
    cell=cell2mat(struct2cell(load(file)));
    tumor(:,i) = cell(:,25)+cell(:,26)+cell(:,27)+cell(:,28);
    tumor(:,i) = ((3*tumor(:,i).*2.572e-9.*0.8)./(0.37*4*3.14)).^(1/3);
    ratio(:,i) = (tumor(:,i)- tumor_0(i))./tumor_0(i)*100;
end

%% Plot

figure(1)
for i = 1:sample
    if ratio(jianc,i) >= 20
        A(1) = A(1)+1;
        plot(1:T:day,ratio(1:T:day,i),'LineWidth',L,'Color',cell2mat(color(3,1)),'DisplayName','PD');
        hold on
    elseif -30 <= min(ratio(jianc:day,i)) && min(ratio(jianc:day,i)) < 20
        A(2) = A(2)+1;
        plot(1:T:day,ratio(1:T:day,i),'LineWidth',L,'Color',cell2mat(color(3,4)),'DisplayName','SD');
        hold on
    elseif  -80 <= min(ratio(jianc:day,i)) &&min(ratio(jianc:day,i)) <-30
        A(3) = A(3)+1;
        plot(1:T:day,ratio(1:T:day,i),'LineWidth',L,'Color',cell2mat(color(3,10)),'DisplayName','PR/CR');
        hold on
    elseif  min(ratio(jianc:day,i)) <-80
        A(4) = A(4)+1;
        plot(1:T:day,ratio(1:T:day,i),'LineWidth',L,'Color',cell2mat(color(3,13)),'DisplayName','PR/CR');
        hold on
    end    
    plot([0,day],[20,20],'--k');
    plot([0,day],[-30,-30],'--k');
end

set(gca,'XLim',[0,day]);
set(gca,'YLim',[-100,100]);
set(gca,'Linewidth',1);   % Set outer border
set(gca,'FontSize',word_size);   % Set font size
set(gca,'FontName','Airal');  % Set font name
xlabel('time (day)');
ylabel('Percentage change in tumor (%)'); 
title('Atezolizumab + Bevacizumab')
set(gcf,'unit','centimeters','position',[1 1 lengthe hight]); 
print('-dpng','-r600','picture\Percentage_change_PCPS.png')


end