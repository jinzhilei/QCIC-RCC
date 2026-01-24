
function Picture_PK_model_atezolizumab_central_peripheral()

sample = 1000;

for i=1:1

    file_sample=['./VP/','VPsample',num2str(i),'.mat'];
    A=cell2mat(struct2cell(load(file_sample)));
    plot(A(1:210,2),A(1:210,1),'Color',[50 205 50]./255);
    hold on

end

cycle1 = 'cycle 1';
text(30,750,cycle1,'FontSize',14,'Interpreter','latex');   
annotation('arrow','X',[0.225,0.275],'Y',[0.56,0.65]);   

cycle2 = 'cycle 2';
text(220,750,cycle2,'FontSize',14,'Interpreter','latex');   
annotation('arrow','X',[0.535,0.5],'Y',[0.56,0.65]); 

cycle3 = 'cycle 3';
text(50,1400,cycle3,'FontSize',14,'Interpreter','latex');   
annotation('arrow','X',[0.275,0.36],'Y',[0.845,0.81]); 

cycle4 = 'limit cycle';
text(280,1250,cycle4,'FontSize',14,'Interpreter','latex');   
annotation('arrow','X',[0.65,0.55],'Y',[0.76,0.69]); 


%% Output image setting

% Set coordinate range
set(gca,'XLim',[0 450]);
set(gca,'YLim',[0 1500]);
% Set outer border
set(gca,'Linewidth',1);
% Set font size
set(gca,'FontSize',14);
% Set font name
set(gca,'FontName','Times New Roman');
% Set label
xlabel('Atezolizumab peripheral concentration ($\mu$g/mL)','Interpreter','latex');
ylabel('Atezolizumab serum concentration ($\mu$g/mL)','Interpreter','latex');
% Set the image name, format, and resolution
set(gcf,'unit','centimeters','position',[10 10 15 12])
print('Picture_Atezolizumab_central_peripheral','-dpng','-r600')