function picture(sample,r)
tumor_0 = ((3*r.*2.572e-9.*0.8)./(0.37*4*3.14)).^(1/3);
for i =1:sample
    file=['./VP_A/','VPsample_A',num2str(i),'.mat'];
    cell=cell2mat(struct2cell(load(file)));
    tumor_1(:,i) = ((3*cell(:,25).*2.572e-9.*0.8)./(0.37*4*3.14)).^(1/3);
    tumor_2(:,i) = ((3*cell(:,26).*2.572e-9.*0.8)./(0.37*4*3.14)).^(1/3);
    ratio_1(:,i) = (tumor_1(:,i)- tumor_0(i))./tumor_0(i)*100;
    ratio_2(:,i) = (tumor_2(:,i)- tumor_0(i))./tumor_0(i)*100;
    drug(:,i) = cell(:,29);
    cell_Ratio_1(:,i) = cell(:,25) ;
    cell_Ratio_2(:,i) = cell(:,26);
end
figure(2)
subplot(2,2,1)
for i = 1:sample
    plot(1:1:360,ratio_1(1:1:360,i));
    hold on
end
plot([43,43],[-100,100],'--k');
set(gca,'XLim',[0,360]);
set(gca,'YLim',[-100,100]);
subplot(2,2,2)
for i = 1:sample
    plot(1:1:360,ratio_2(1:1:360,i));
    hold on
end
plot([43,43],[-100,100],'--k');
set(gca,'XLim',[0,360]);
set(gca,'YLim',[-100,100]);


end