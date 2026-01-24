function picture(sample)
day = 360;
file_R='./VP_AB/VP_tumor_number.mat';
r=cell2mat(struct2cell(load(file_R)));
tumor_0 = ((3*r.*2.572e-9.*0.8)./(0.37*4*3.14)).^(1/3);
T=21;
for i =1:sample
    file=['./VP_AB/','VPsample_AB',num2str(i),'.mat'];
    cell=cell2mat(struct2cell(load(file)));
    tumor_1(:,i) = ((3*cell(:,25).*2.572e-9.*0.8)./(0.37*4*3.14)).^(1/3);
    tumor_2(:,i) = ((3*cell(:,26).*2.572e-9.*0.8)./(0.37*4*3.14)).^(1/3);
    tumor_3(:,i) = ((3*cell(:,27).*2.572e-9.*0.8)./(0.37*4*3.14)).^(1/3);
    tumor_4(:,i) = ((3*cell(:,28).*2.572e-9.*0.8)./(0.37*4*3.14)).^(1/3);
    ratio_1(:,i) = (tumor_1(:,i)- tumor_0(i))./tumor_0(i)*100;
    ratio_2(:,i) = (tumor_2(:,i)- tumor_0(i))./tumor_0(i)*100;
    ratio_3(:,i) = (tumor_3(:,i)- tumor_0(i))./tumor_0(i)*100;
    ratio_4(:,i) = (tumor_4(:,i)- tumor_0(i))./tumor_0(i)*100;
    drug_A(:,i) = cell(:,29);
    drug_B(:,i) = cell(:,31);
end
figure(2)
subplot(2,2,1)
for i = 1:sample
    plot(1:T:360,ratio_1(1:T:360,i));
    hold on
    plot([0,day],[20,20],'--k');
    plot([0,day],[-30,-30],'--k');
    plot([0,day],[-80,-80],'--k');
end
set(gca,'XLim',[0,360]);
set(gca,'YLim',[-100,100]);
xlabel('X1');

subplot(2,2,2)
for i = 1:sample
    plot(1:T:360,ratio_2(1:T:360,i));
    hold on
    plot([0,day],[20,20],'--k');
    plot([0,day],[-30,-30],'--k');
    plot([0,day],[-80,-80],'--k');
end
set(gca,'XLim',[0,360]);
set(gca,'YLim',[-100,100]);
xlabel('X2');

subplot(2,2,3)
for i = 1:sample
    plot(1:T:360,ratio_3(1:T:360,i));
    hold on
    plot([0,day],[20,20],'--k');
    plot([0,day],[-30,-30],'--k');
    plot([0,day],[-80,-80],'--k');
end
set(gca,'XLim',[0,360]);
set(gca,'YLim',[-100,100]);
xlabel('X3');

subplot(2,2,4)
for i = 1:sample
    plot(1:T:360,ratio_4(1:T:360,i));
    hold on
    plot([0,day],[20,20],'--k');
    plot([0,day],[-30,-30],'--k');
    plot([0,day],[-80,-80],'--k');
end
set(gca,'XLim',[0,360]);
set(gca,'YLim',[-100,100]);
xlabel('X4');





end