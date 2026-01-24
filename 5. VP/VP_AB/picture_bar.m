function picture_bar()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This module is used to draw waterfall diagrams
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sample= 100;
jiance = 43;
day = 360;
lengthe = 8*1.8;
hight = 6;
word_size = 9;
color = ColorMatrix();
file_R='./VP_AB/VP_tumor_number.mat';
r=cell2mat(struct2cell(load(file_R)));
tumor = zeros(day,sample);
ratio_all = zeros(day,sample);
ratio = zeros(1,sample);
tumor_0 = ((r.*2.572e-9.*0.8*3)./(0.37*4*3.14)).^(1/3);

%% Calculation of tumor evaluation indicators

for i =1:sample
    file=['./VP_AB/','VPsample_AB',num2str(i),'.mat'];
    cell=cell2mat(struct2cell(load(file)));
    tumor(:,i) = cell(:,25)+cell(:,26)+cell(:,27)+cell(:,28);
    tumor(:,i) = ((3*tumor(:,i).*2.572e-9.*0.8)./(0.37*4*3.14)).^(1/3);
    ratio_all(:,i) = (tumor(:,i)- tumor_0(i))./tumor_0(i)*100;
    ratio(i) =  min(ratio_all(jiance:day,i));
end

[~, I_sort] = sort(ratio,'descend');
S_sorted = ratio(I_sort);

%% Plot

figure(2)
b=bar(S_sorted,0.7);
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
plot([21.5,21.5],[-100,100],'--k');
set(gca,'XLim',[0,sample+9]);
set(gca,'YLim',[-100,100]);
set(gca,'Linewidth',1);   % Set outer border
set(gca,'FontSize',word_size);   % Set font size
set(gca,'FontName','Airal');  % Set font name
ylabel('Percentage change in tumor (%)');
title('Atezolizumab + Bevacizumab')
text(103,60,{'PD'},'FontSize',word_size);
text(103,-5,{'SD'},'FontSize',word_size);
text(101,-50,{'PR/CR'},'FontSize',word_size);
set(gcf,'unit','centimeters','position',[1 1 lengthe hight]); 
print('-dpng','-r600','picture\Percentage_change_in_tumor_size.png')
end