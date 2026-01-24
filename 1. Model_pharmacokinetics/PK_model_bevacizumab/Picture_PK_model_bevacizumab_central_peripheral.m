
function Picture_PK_model_bevacizumab_central_peripheral()

sample = 1000;

for i=1:1

    file_sample=['./VP/','VPsample',num2str(i),'.mat'];
    A=cell2mat(struct2cell(load(file_sample)));
    plot(A(1:210,2),A(1:210,1),'Color',[30 144 255]./255);
    hold on

end

cycle1 = 'cycle 1';
text(20,740,cycle1,'FontSize',14,'Interpreter','latex');   
annotation('arrow','X',[0.2,0.25],'Y',[0.6,0.68]);   

cycle2 = 'cycle 2';
text(200,750,cycle2,'FontSize',14,'Interpreter','latex');   
annotation('arrow','X',[0.5,0.45],'Y',[0.6,0.675]); 

cycle3 = 'cycle 3';
text(20,1220,cycle3,'FontSize',14,'Interpreter','latex');   
annotation('arrow','X',[0.225,0.295],'Y',[0.85,0.81]); 

cycle4 = 'limit cycle';
text(280,1000,cycle4,'FontSize',14,'Interpreter','latex');   
annotation('arrow','X',[0.67,0.6],'Y',[0.7,0.626]); 


%% Output image setting

% Set coordinate range
set(gca,'XLim',[0 450]);
set(gca,'YLim',[0 1300]);
% Set outer border
set(gca,'Linewidth',1);
% Set font size
set(gca,'FontSize',14);
% Set font name
set(gca,'FontName','Times New Roman');
% Set label
xlabel('Bevacizumab peripheral concentration ($\mu$g/mL)','Interpreter','latex');
ylabel('Bevacizumab serum concentration ($\mu$g/mL)','Interpreter','latex');
% Set the image name, format, and resolution
set(gcf,'unit','centimeters','position',[10 10 15 12])
print('Picture_Bevacizumab_central_peripheral','-dpng','-r600')